// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2014, 2026 Travis Stenborg

#include <mpi.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem> // NOLINT.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

using json = nlohmann::json;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

// Include "iomanip" if setprecision calls are enabled.
// #include <iomanip>

/*
A "Stable ID" is a uint64_t variable used to identify covers.
It's a platform-independent 64-bit (8 byte) integer ID, i.e.:

bits   00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000

bit #  6      5        4        4        3        2        1        8        0
       3      6        8        0        2        4        6


To build deterministic IDs for covers, this bit layout adopted is:

  Singleton namespace : bit 63   A flag specifying a singleton cover.
  Target One :    bits 32-62     Accommodates 2^31 = 2,147,483,648 max targets.
  Target Two :    bits  1-31     Accommodates 2^31 = 2,147,483,648 max targets.
  Cover Variant : bit 0          1 or 2 variants are used. 2^1 = 2 values is ok.

The IDs are independent of MPI rank and execution order, allowing deterministic
tie-breaking when multiple covers contain the same number of targets, hence a
"Stable ID".

N.B. The Stable ID encoding supports far more targets than the current
     implementation, which stores target indices as int32_t.
*/

constexpr int STABLE_ID_BIT_SHIFT_1 = 1;
constexpr int STABLE_ID_BIT_SHIFT_32 = 32;

struct CandidateCover {
  // Store positions in degrees.
  double GalacticLat = 0.0, GalacticLong = 0.0;
  // Platform-independent 64-bit (8 byte) integer ID.
  uint64_t StableID = 0;
  int TargetsCovered = 0;
};

struct GreedyCover {
  uint64_t StableID = UINT64_MAX;
  int TargetsCovered = 0;
};

struct IndexPair {
  // Store separations in degrees.
  double SeparationLat = 0.0, SeparationLong = 0.0;
  int32_t TargetOne = 0, TargetTwo = 0;
};

struct Target {
  // Store positions in degrees.
  double GalacticLat = 0.0, GalacticLong = 0.0;
  bool Covered = false;
};

void GetWorkRange(int intTotalWork, int intRank, int intSize, int &intStart,
                  int &intEnd) {
  // A helper function used to help distribute work evenly across processes.
  // It aims to ensure process of rank 0 always has work, if any exists.

  // "/" is integer division; divide and discard the remainder. E.g., 5/2 = 2.
  // "%" is the modulus operator. It gets the remainder. E.g., 5%2 = 1.
  int intBaseWork = intTotalWork / intSize;
  int intRemainder = intTotalWork % intSize;

  if (intRank < intRemainder) {
    intStart = intRank * (intBaseWork + 1);
    intEnd = intStart + (intBaseWork + 1);
  } else {
    intStart = intRemainder * (intBaseWork + 1) +
               (intRank - intRemainder) * intBaseWork;
    intEnd = intStart + intBaseWork;
  }
}

vector<Target> ReadTargets(const std::string &strFileName) {
  // This function reads target coordinate data from a JSON file, and returns
  // it in a vector of Target structs.
  //
  // Assume the JSON data are formatted as follows:
  //
  // {
  //   "galactic_latitude_longitude": [
  //     [-5.03890867670439, -2.30821603331922],
  //     [4.99738096958076, 3.96327582850898]
  //   ]
  // }

  // Create an input file stream "file" to read a specified file.
  std::ifstream file(strFileName);
  if (!file) {
    throw std::runtime_error("Could not open JSON file.");
  }

  // Read JSON data from the specified file.
  const json data = json::parse(file);
  const auto &coordinates = data.at("galactic_latitude_longitude");

  // Create a vector to hold the targets, and reserve memory to hold them.
  vector<Target> vtTargets;
  vtTargets.reserve(coordinates.size());

  for (const auto &pair : coordinates) {
    if (!pair.is_array() || pair.size() != 2) {
      throw std::runtime_error(
          "Each coordinate must be an array of Galactic [latitude, longitude]");
    }

    Target target;
    target.GalacticLat = pair.at(0).get<double>();
    target.GalacticLong = pair.at(1).get<double>();
    vtTargets.push_back(target);
  }

  return vtTargets;
}

bool IsEpsilonClose(double a, double b) {
  // This function tests if two input values are equal, to within the bounds of
  // a specific epsilon.

  // The data used here have a 14 decimal place precision.
  // The data range is -5.03890867670439 to 4.99738096958076.
  // It's assumed an epsilon of 5e-15 is suitable for comparing such data.
  constexpr double EPSILON = 5e-15;

  return std::abs(a - b) <= EPSILON;
}

int main(int argc, char *argv[]) {
  constexpr double FOV_RADIUS = 1.05;
  constexpr double FOV_RADIUS_POW2 = FOV_RADIUS * FOV_RADIUS;
  constexpr double FOV_DIAMETER = 2 * FOV_RADIUS;

  vector<Target> vtTarget;
  vector<IndexPair> vtIndexPair;

  vector<CandidateCover> vtCandidateCover;
  vector<CandidateCover> vtFilteredCandidateCover;
  vector<CandidateCover> vtGlobalCandidateCover;
  CandidateCover CandidateCoverExtentProbe[2];
  CandidateCover CandidateCoverFiller;

  GreedyCover GreedyCoverGlobal;
  GreedyCover GreedyCoverLocal;
  vector<GreedyCover> vtGreedyCovers;

  // Initialise MPI datatypes and related management variables.
  // N.B. MPI_Aint is an integer data type used to hold a memory address or
  //      byte displacement.
  MPI_Datatype typesArray[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_UINT64_T};
  MPI_Aint offsetsArray[4];
  int intBlocklengthsArray[4] = {1, 1, 1, 1};
  MPI_Datatype MPI_CandidateCover, MPI_CandidateCoverPrototype;
  MPI_Aint a0, a1, offsetLowerBound, offsetExtent;

  std::filesystem::path pthDataFile;
  std::string strDataFile;
  double dblApothem = 0.0;
  double dblAvgLat = 0.0;
  double dblAvgLong = 0.0;
  double dblScaledSeparation = 0.0;
  double dblScaledUnitPerpVectorLat = 0.0;
  double dblScaledUnitPerpVectorLong = 0.0;
  double dblSeparation = 0.0;
  double dblSeparationLat = 0.0;
  double dblSeparationLong = 0.0;
  double dblTimeCompleted = 0.0;
  double dblTimeStarted = 0.0;
  uint64_t ui64CurrentID = 0;
  uint64_t ui64GlobalStableID = 0;
  uint64_t ui64LocalBestStableID = 0;
  uint64_t ui64StableID1 = 0;
  uint64_t ui64StableID2 = 0;
  std::size_t uiNumPairs = 0.0;
  int32_t i32Covers = 0;
  int32_t i32MaxCovers = 0;
  int32_t i32TotalCovers = 0;
  int intEnd = 0;
  int intGlobalMaxCoverIndex = 0;
  int intI = 0;
  int intItems = 0;
  int intJ = 0;
  int intMaxCoverIndex = 0;
  int intNumTargets = 0;
  int intRank = 0;
  int intSize = 0;
  int intStart = 0;
  int intTargetsMax = 0;
  int intTargetsTemp = 0;
  int intX = 0;

  // Ensure an input data file to process has been supplied.
  if (argc > 1) {
    strDataFile = argv[1];
  } else {
    cerr << "Usage:\n" << argv[0] << " <input-file.json>\n";
    return 1;
  }

  // Ensure the specified input data file actually exists.
  pthDataFile = strDataFile;
  if (!std::filesystem::is_regular_file(pthDataFile)) {
    cerr << "Input file not found or is not a regular file: " << pthDataFile
         << "\n";
    return 1;
  }

  // Initialise target vector and supporting variables.
  vtTarget = ReadTargets(strDataFile);
  intNumTargets = vtTarget.size();
  uiNumPairs =
      static_cast<std::size_t>(intNumTargets) * (intNumTargets - 1) / 2;
  // Set the number of index pairs to 0.
  vtIndexPair.clear();
  // Reserve memory for uiNumPairs index pairs, without initialising any.
  vtIndexPair.reserve(uiNumPairs);

  // Initialise variables.
  CandidateCoverFiller.GalacticLat = 0.0;
  CandidateCoverFiller.GalacticLong = 0.0;
  CandidateCoverFiller.TargetsCovered = -1;

  // Initialise MPI.
  MPI_Init(&argc, &argv);

  // Get the number of processes available in the communicator.
  MPI_Comm_size(MPI_COMM_WORLD, &intSize);
  // Get a process's rank within the communicator, from 0 to (intSize - 1).
  MPI_Comm_rank(MPI_COMM_WORLD, &intRank);

  dblTimeStarted = MPI_Wtime();

  // Create an MPI struct for the CandidateCover struct.
  intItems = std::size(intBlocklengthsArray);
  offsetsArray[0] = offsetof(CandidateCover, GalacticLat);
  offsetsArray[1] = offsetof(CandidateCover, GalacticLong);
  offsetsArray[2] = offsetof(CandidateCover, TargetsCovered);
  offsetsArray[3] = offsetof(CandidateCover, StableID);
  MPI_Type_create_struct(intItems, intBlocklengthsArray, offsetsArray,
                         typesArray, &MPI_CandidateCoverPrototype);

  // Resize the MPI struct.
  // Ensure its length matches its length implemented in memory.
  MPI_Type_get_extent(MPI_CandidateCoverPrototype, &offsetLowerBound,
                      &offsetExtent);

  // Calculate the vector's inter-element distance.
  MPI_Get_address(&CandidateCoverExtentProbe[0], &a0);
  MPI_Get_address(&CandidateCoverExtentProbe[1], &a1);
  offsetExtent = MPI_Aint_diff(a1, a0);

  // MPI_Type_create_struct describes the members in a struct, but not the
  // spacing between consecutive vector elements. Therefore, use
  // MPI_Type_create_resized to create a resized datatype whose extent matches
  // the actual in-memory distance between adjacent CandidateCover objects,
  // ensuring MPI communication is portable across compilers and architectures.
  MPI_Type_create_resized(MPI_CandidateCoverPrototype, offsetLowerBound,
                          offsetExtent, &MPI_CandidateCover);
  MPI_Type_commit(&MPI_CandidateCover);

  do {
    // Initialise a target pair indices array.
    for (intI = 0; intI <= (intNumTargets - 2); intI++) {
      for (intJ = (intI + 1); intJ <= (intNumTargets - 1); intJ++) {
        // Append an element to vector vtIndexPair.
        // N.B. An IndexPair struct has the following member variables:
        //          SeparationLat, SeparationLong, TargetOne, TargetTwo.
        vtIndexPair.push_back({0.0, 0.0, intI, intJ});
      }
    }

    // Configure parallelisation of target pair proximity comparisons and cover
    // construction.
    // N.B. Configure a unique "intStart" and "intEnd" for each rank.
    GetWorkRange(static_cast<int>(vtIndexPair.size()), intRank, intSize,
                 intStart, intEnd);

    // Calculate target separations and assign covers.
    // Note, workload here is parallelised.
    i32Covers = 0;
    vtCandidateCover.clear();
    for (intX = intStart; intX < intEnd; intX++) {
      // Some separations may be stored as negative values.
      // This is unphysical, but they are soon squared and never reused.
      vtIndexPair.at(intX).SeparationLat =
          vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLat -
          vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLat;
      vtIndexPair.at(intX).SeparationLong =
          vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLong -
          vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLong;

      dblSeparation = sqrt(vtIndexPair.at(intX).SeparationLat *
                               vtIndexPair.at(intX).SeparationLat +
                           vtIndexPair.at(intX).SeparationLong *
                               vtIndexPair.at(intX).SeparationLong);

      if (IsEpsilonClose(dblSeparation, FOV_DIAMETER)) {
        // Add one cover for two targets.
        vtCandidateCover.resize(i32Covers + 1);

        // Calculate the midpoint of the chord between the targets.
        vtCandidateCover.at(i32Covers).GalacticLat =
            0.5 * (vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLat +
                   vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLat);
        vtCandidateCover.at(i32Covers).GalacticLong =
            0.5 * (vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLong +
                   vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLong);
        // Update the cover ID.
        // N.B. "<<" is a bitwise left shift.
        vtCandidateCover.at(i32Covers).StableID =
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetOne)
             << STABLE_ID_BIT_SHIFT_32) |
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetTwo)
             << STABLE_ID_BIT_SHIFT_1);

        i32Covers += 1;
      } else if (dblSeparation < FOV_DIAMETER) {
        // Add two covers.

        // Calculate the apothem; the distance from a cover centre to the
        // midpoint of the chord between the targets.
        dblApothem =
            sqrt(FOV_RADIUS_POW2 - (0.25 * dblSeparation * dblSeparation));

        // Calculate the midpoint of the chord between the targets.
        dblAvgLat =
            0.5 * (vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLat +
                   vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLat);
        dblAvgLong =
            0.5 * (vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLong +
                   vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLong);

        // Calculate distances from the chord midpoint to the cover centre.
        dblScaledSeparation = dblApothem / dblSeparation;
        dblScaledUnitPerpVectorLat =
            dblScaledSeparation * vtIndexPair.at(intX).SeparationLat;
        dblScaledUnitPerpVectorLong =
            dblScaledSeparation * vtIndexPair.at(intX).SeparationLong;

        // Cover 1.
        // N.B. Keep the comments in this section, or autoformatters mangle it.
        vtCandidateCover.resize(i32Covers + 1);
        // Update latitude.
        vtCandidateCover.at(i32Covers).GalacticLat =
            dblAvgLat + dblScaledUnitPerpVectorLong;
        // Update longitude.
        vtCandidateCover.at(i32Covers).GalacticLong =
            dblAvgLong - dblScaledUnitPerpVectorLat;
        // Update the cover ID.
        ui64StableID1 = (static_cast<uint64_t>(vtIndexPair.at(intX).TargetOne)
                         << STABLE_ID_BIT_SHIFT_32) |
                        (static_cast<uint64_t>(vtIndexPair.at(intX).TargetTwo)
                         << STABLE_ID_BIT_SHIFT_1) |
                        0;
        vtCandidateCover.at(i32Covers).StableID = ui64StableID1;
        // Increase cover counter.
        i32Covers += 1;

        // Cover 2.
        // N.B. Keep the comments in this section, or autoformatters mangle it.
        vtCandidateCover.resize(i32Covers + 1);
        // Update latitude.
        vtCandidateCover.at(i32Covers).GalacticLat =
            dblAvgLat - dblScaledUnitPerpVectorLong;
        // Update longitude.
        vtCandidateCover.at(i32Covers).GalacticLong =
            dblAvgLong + dblScaledUnitPerpVectorLat;
        // Update the cover ID.
        ui64StableID2 = (static_cast<uint64_t>(vtIndexPair.at(intX).TargetOne)
                         << STABLE_ID_BIT_SHIFT_32) |
                        (static_cast<uint64_t>(vtIndexPair.at(intX).TargetTwo)
                         << STABLE_ID_BIT_SHIFT_1) |
                        1;
        vtCandidateCover.at(i32Covers).StableID = ui64StableID2;
        // Increase cover counter.
        i32Covers += 1;
      }
    }

    // Cover singleton targets.
    for (intX = 0; intX < intNumTargets; intX++) {
      vtCandidateCover.resize(i32Covers + 1);
      // Update latitude.
      vtCandidateCover.at(i32Covers).GalacticLat =
          vtTarget.at(intX).GalacticLat;
      // Update longitude.
      vtCandidateCover.at(i32Covers).GalacticLong =
          vtTarget.at(intX).GalacticLong;

      // Update the 64-bit cover ID.
      // N.B. Singleton covers IDs go in a separate namespace by setting the
      // most-significant bit, to prevent collision with pair-derived cover IDs.
      // N.B. 1ULL is the value 1 stored as an unsigned long long.
      //      1ULL << 63 = ...
      // 10000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
      vtCandidateCover.at(i32Covers).StableID =
          (1ULL << 63) | static_cast<uint64_t>(intX);
      // Increase cover counter.
      i32Covers += 1;
    }

    // Free local index pair data.
    // N.B. clear() removes all vector elements, leaving it with a size of 0.
    vtIndexPair.clear();

    // Determine the local maximum number of candidate covers created.
    // N.B. MPI_Allreduce is a synchronising collective communication function.
    // It doesn't need leading or trailing barriers.
    MPI_Allreduce(&i32Covers, &i32MaxCovers, 1, MPI_INT32_T, MPI_MAX,
                  MPI_COMM_WORLD);
    if (i32MaxCovers == 0) {
      if (intRank == 0) {
        cout << "Remaining targets per pointing <= 1." << endl;
      }
      break;
    }

    // Give each candidate cover a uniform extent to facilitate Allgather.
    if (static_cast<std::size_t>(vtCandidateCover.size()) < i32MaxCovers) {
      // Set the number of elements to i32MaxCovers.
      // Initialise any new elements with CandidateCoverFiller.
      vtCandidateCover.resize(i32MaxCovers, CandidateCoverFiller);
    }

    // Build vector to receive global candidate cover data into.
    vtGlobalCandidateCover.resize(i32MaxCovers * intSize);

    // Distribute local candidate cover data across all ranks.
    // N.B. ".data()" returns a pointer to the underlying raw array used by the
    // vectors to store their data.
    // N.B. MPI_Allgather is a synchronising collective communication function.
    // It doesn't need leading or trailing barriers.
    MPI_Allgather(vtCandidateCover.data(), i32MaxCovers, MPI_CandidateCover,
                  vtGlobalCandidateCover.data(), i32MaxCovers,
                  MPI_CandidateCover, MPI_COMM_WORLD);

    // Free unused candidate covers data.
    // N.B. clear() removes all vector elements, leaving it at size 0.
    vtCandidateCover.clear();

    // Reserve memory for filtered candidate covers, without initialising any.
    vtFilteredCandidateCover.clear();
    vtFilteredCandidateCover.reserve(vtGlobalCandidateCover.size());

    // Remove any default filler candidate covers.
    for (const auto &cover : vtGlobalCandidateCover) {
      if (cover.TargetsCovered != -1)
        // Add "cover" to vtFilteredCandidateCover.
        vtFilteredCandidateCover.push_back(cover);
    }
    vtGlobalCandidateCover.swap(vtFilteredCandidateCover);
    i32TotalCovers = vtGlobalCandidateCover.size();

    // Configure parallelisation of target coverage calculations.
    GetWorkRange(i32TotalCovers, intRank, intSize, intStart, intEnd);

    // Prepare for deterministic greedy set cover tie-breaking.
    // A "stable ID" is used for consistent cover selection, regardless of
    // rank count, workload partitioning, or candidate cover ordering.
    ui64LocalBestStableID = UINT64_MAX;
    intMaxCoverIndex = -1;

    // Reset cover counts.
    for (intI = 0; intI < i32TotalCovers; intI++) {
      vtGlobalCandidateCover.at(intI).TargetsCovered = 0;
    }

    // Get each rank to determine the greediest cover within its assigned
    // subset of the candidate covers.
    // Note, workload here is parallelised.
    intTargetsMax = 0;
    for (intI = intStart; intI < intEnd; intI++) {
      intTargetsTemp = 0;

      for (intJ = 0; intJ < intNumTargets; intJ++) {
        // Some separations may be stored as negative values.
        // This is unphysical, but they are soon squared and never reused.
        dblSeparationLat = vtTarget.at(intJ).GalacticLat -
                           vtGlobalCandidateCover.at(intI).GalacticLat;
        dblSeparationLong = vtTarget.at(intJ).GalacticLong -
                            vtGlobalCandidateCover.at(intI).GalacticLong;

        // A target is covered if it's within FOV_RADIUS of the cover's centre.
        //
        // N.B. Use squared-distance comparison for performance:
        // Slower: if (sqrt(a*a + b*b) < FOV_RADIUS)
        // Faster: if ((a*a + b*b) < FOV_RADIUS_POW2)
        //
        if ((dblSeparationLat * dblSeparationLat +
             dblSeparationLong * dblSeparationLong) < FOV_RADIUS_POW2) {
          vtGlobalCandidateCover.at(intI).TargetsCovered += 1;
          intTargetsTemp += 1;
        }
      }
      // Assess greediest covers, with robust tie-breaking.
      ui64CurrentID = vtGlobalCandidateCover.at(intI).StableID;
      if (intTargetsTemp > intTargetsMax ||
          (intTargetsTemp == intTargetsMax &&
           ui64CurrentID < ui64LocalBestStableID)) {
        intTargetsMax = intTargetsTemp;
        intMaxCoverIndex = intI;
        ui64LocalBestStableID = ui64CurrentID;
      }
    }

    // Collect each rank's GreedyCoverLocal into vtGreedyCovers.
    // N.B. MPI_Allgather is a synchronising collective communication function.
    // It doesn't need leading or trailing barriers.
    GreedyCoverLocal.TargetsCovered = intTargetsMax;
    GreedyCoverLocal.StableID = ui64LocalBestStableID;
    vtGreedyCovers.resize(intSize);
    MPI_Allgather(&GreedyCoverLocal, sizeof(GreedyCover), MPI_BYTE,
                  vtGreedyCovers.data(), sizeof(GreedyCover), MPI_BYTE,
                  MPI_COMM_WORLD);

    GreedyCoverGlobal.TargetsCovered = 0;
    GreedyCoverGlobal.StableID = UINT64_MAX;

    // Determine the globally greediest cover.
    //
    // Primary criterion:   maximum targets covered.
    // Secondary criterion: smallest StableID.
    //
    // A StableID tie-break guarantees all MPI configurations produce identical
    // cover selections when multiple covers enclose the same number of targets.
    for (const auto &GreedyCoverCandidate : vtGreedyCovers) {
      if (GreedyCoverCandidate.TargetsCovered >
              GreedyCoverGlobal.TargetsCovered ||
          (GreedyCoverCandidate.TargetsCovered ==
               GreedyCoverGlobal.TargetsCovered &&
           GreedyCoverCandidate.StableID < GreedyCoverGlobal.StableID)) {
        GreedyCoverGlobal = GreedyCoverCandidate;
      }
    }

    // Locate the greediest cover within the global candidate cover list.
    // Stable IDs are unique, so can be used for an exact match.
    ui64GlobalStableID = GreedyCoverGlobal.StableID;
    intGlobalMaxCoverIndex = -1;
    for (intI = 0; intI < i32TotalCovers; intI++) {
      if (vtGlobalCandidateCover.at(intI).StableID == ui64GlobalStableID) {
        intGlobalMaxCoverIndex = intI;
        break;
      }
    }
    // Safety check.
    if (intGlobalMaxCoverIndex < 0) {
      cerr << "Error: Failed to locate global greediest cover." << endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    intTargetsMax = GreedyCoverGlobal.TargetsCovered;

    // Remove covered targets.
    for (intJ = 0; intJ < intNumTargets; intJ++) {
      // Some separations may be stored as negative values.
      // This is unphysical, but they are soon squared and never reused.
      dblSeparationLat =
          vtTarget.at(intJ).GalacticLat -
          vtGlobalCandidateCover.at(intGlobalMaxCoverIndex).GalacticLat;
      dblSeparationLong =
          vtTarget.at(intJ).GalacticLong -
          vtGlobalCandidateCover.at(intGlobalMaxCoverIndex).GalacticLong;

      // A target is covered if it's within FOV_RADIUS of the cover's centre.
      // N.B. Use squared-distance comparison for performance:
      // Slower: if (sqrt(a*a + b*b) < FOV_RADIUS)
      // Faster: if ((a*a + b*b) < FOV_RADIUS_POW2)
      if ((dblSeparationLat * dblSeparationLat +
           dblSeparationLong * dblSeparationLong) < FOV_RADIUS_POW2) {
        vtTarget.at(intJ).Covered = true;
      }
    }

    // Iterate backwards and erase elements from vtTarget.
    // Removing higher-index elements first avoids removing valid elements.
    for (intI = (intNumTargets - 1); intI >= 0; intI--) {
      if (vtTarget.at(intI).Covered == true) {
        vtTarget.erase(vtTarget.begin() + intI);
      }
    }
    intNumTargets = static_cast<int>(vtTarget.size());

    // Display results.
    if (intRank == 0) {
      cout << "Greediest Cover" << endl;
      cout << "   Galactic latitude: "
           << vtGlobalCandidateCover.at(intGlobalMaxCoverIndex).GalacticLat
           << endl;
      cout << "   Galactic longitude: "
           << vtGlobalCandidateCover.at(intGlobalMaxCoverIndex).GalacticLong
           << endl;
      cout << "   Targets covered: " << intTargetsMax << endl;
    }

    // Free vector data.
    vtGlobalCandidateCover.clear();
  } while (intNumTargets > 0);

  MPI_Barrier(MPI_COMM_WORLD);
  dblTimeCompleted = MPI_Wtime();

  // Terminate the MPI execution environment.
  MPI_Finalize();

  cout << "MPI Execution time: " << dblTimeCompleted - dblTimeStarted << endl;

  return 0;
}

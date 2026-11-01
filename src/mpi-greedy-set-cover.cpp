// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2014, 2026 Travis Stenborg

#include <mpi.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

// Include "iomanip" if setprecision calls are enabled.
// #include <iomanip>

struct IndexPair {
  // Store separations in degrees.
  double SeparationLat = 0.0, SeparationLong = 0.0;
  int16_t TargetOne = 0, TargetTwo = 0;
};

struct CoveringCircle {
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

bool isEpsilonClose(double a, double b) {
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
  int intNumTargets = 1100;

  const std::size_t numPairs =
      static_cast<std::size_t>(intNumTargets) * (intNumTargets - 1) / 2;

  vector<Target> vtTarget(intNumTargets);
  vector<IndexPair> vtIndexPair;
  vtIndexPair.clear();
  vtIndexPair.reserve(numPairs);

  vector<CoveringCircle> vtCoveringCircle;
  vector<CoveringCircle> vtGlobalCoveringCircle;
  CoveringCircle CoveringCircleExtentProbe[2];
  CoveringCircle CoveringCircleFiller;

  GreedyCover GreedyCoverGlobal;
  GreedyCover GreedyCoverLocal;
  vector<GreedyCover> vtGreedyCovers;

  // Segmentation fault without some of these initialisations.
  // N.B. MPI_Aint is an integer data type used to hold a memory address or
  //      byte displacement.
  MPI_Datatype typesArray[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_UINT64_T};
  MPI_Aint offsetsArray[4];
  int intBlocklengthsArray[4] = {1, 1, 1, 1};
  MPI_Datatype MPI_CoveringCircle, MPI_CoveringCirclePrototype;
  MPI_Aint a0, a1, offsetLowerBound, offsetExtent;

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
  uint64_t ui64StableID12 = 0;
  int32_t i32Covers = 0;
  int32_t i32MaxCovers = 0;
  int32_t i32TotalCovers = 0;
  int intEnd = 0;
  int intGlobalMaxCoverIndex = 0;
  int intI = 0;
  int intItems = 0;
  int intJ = 0;
  int intMaxCoverIndex = 0;
  int intRank = 0;
  int intSize = 0;
  int intStart = 0;
  int intTargetsMax = 0;
  int intTargetsTemp = 0;
  int intX = 0;

  // Initialise target array.
  vtTarget.at(0).GalacticLat = -3.43527864026376;
  vtTarget.at(1).GalacticLat = -3.42931387831587;
  vtTarget.at(2).GalacticLat = -3.77850429184769;
  vtTarget.at(3).GalacticLat = -3.41294869905441;
  vtTarget.at(4).GalacticLat = -3.48142759473224;
  vtTarget.at(5).GalacticLat = -4.3025065754045;
  vtTarget.at(6).GalacticLat = -4.33665788789297;
  vtTarget.at(7).GalacticLat = -4.38535270402158;
  vtTarget.at(8).GalacticLat = -3.87704808895843;
  vtTarget.at(9).GalacticLat = -4.55669244846708;
  vtTarget.at(10).GalacticLat = -4.59123555143811;
  vtTarget.at(11).GalacticLat = -4.63509436985048;
  vtTarget.at(12).GalacticLat = -4.62632991346207;
  vtTarget.at(13).GalacticLat = -4.83095769041197;
  vtTarget.at(14).GalacticLat = -3.4531833521599;
  vtTarget.at(15).GalacticLat = -3.5796283725168;
  vtTarget.at(16).GalacticLat = -3.88508484874186;
  vtTarget.at(17).GalacticLat = -4.26111847881099;
  vtTarget.at(18).GalacticLat = -4.61940175592868;
  vtTarget.at(19).GalacticLat = -4.76854528678105;
  vtTarget.at(20).GalacticLat = -1.19459767336516;
  vtTarget.at(21).GalacticLat = -1.20442300353696;
  vtTarget.at(22).GalacticLat = -1.35920351747864;
  vtTarget.at(23).GalacticLat = -1.5985338021377;
  vtTarget.at(24).GalacticLat = -1.81322373104257;
  vtTarget.at(25).GalacticLat = -2.53764997453987;
  vtTarget.at(26).GalacticLat = -2.19150405605689;
  vtTarget.at(27).GalacticLat = -2.7782851374562;
  vtTarget.at(28).GalacticLat = -1.58816032344855;
  vtTarget.at(29).GalacticLat = -1.96853986624778;
  vtTarget.at(30).GalacticLat = -2.02662797953438;
  vtTarget.at(31).GalacticLat = -1.36471246022341;
  vtTarget.at(32).GalacticLat = -1.98459633945926;
  vtTarget.at(33).GalacticLat = -1.80954106739064;
  vtTarget.at(34).GalacticLat = -1.84261027182213;
  vtTarget.at(35).GalacticLat = -2.06712523961646;
  vtTarget.at(36).GalacticLat = -3.21889684661984;
  vtTarget.at(37).GalacticLat = 3.01214960618047;
  vtTarget.at(38).GalacticLat = 2.69086356104855;
  vtTarget.at(39).GalacticLat = 2.24212710825483;
  vtTarget.at(40).GalacticLat = 2.13003004147388;
  vtTarget.at(41).GalacticLat = 2.74790998305413;
  vtTarget.at(42).GalacticLat = 1.95822456165249;
  vtTarget.at(43).GalacticLat = 1.96241778633826;
  vtTarget.at(44).GalacticLat = 1.73191807956738;
  vtTarget.at(45).GalacticLat = 2.00280667072857;
  vtTarget.at(46).GalacticLat = 1.7026908032789;
  vtTarget.at(47).GalacticLat = 2.90174126668037;
  vtTarget.at(48).GalacticLat = 2.06510887870492;
  vtTarget.at(49).GalacticLat = 3.86559554447569;
  vtTarget.at(50).GalacticLat = 3.80381682982785;
  vtTarget.at(51).GalacticLat = 3.71973747508467;
  vtTarget.at(52).GalacticLat = 3.66221220560482;
  vtTarget.at(53).GalacticLat = 2.44862514646353;
  vtTarget.at(54).GalacticLat = 2.36724508243687;
  vtTarget.at(55).GalacticLat = 2.37772883013845;
  vtTarget.at(56).GalacticLat = 2.135901942587;
  vtTarget.at(57).GalacticLat = 3.78985186306797;
  vtTarget.at(58).GalacticLat = 2.86086043247676;
  vtTarget.at(59).GalacticLat = 2.88331353126084;
  vtTarget.at(60).GalacticLat = 2.44058027996283;
  vtTarget.at(61).GalacticLat = 4.18166776173201;
  vtTarget.at(62).GalacticLat = 3.64760816139215;
  vtTarget.at(63).GalacticLat = 4.23263020993331;
  vtTarget.at(64).GalacticLat = 3.74300073376487;
  vtTarget.at(65).GalacticLat = 3.71736089894544;
  vtTarget.at(66).GalacticLat = 3.5894501423651;
  vtTarget.at(67).GalacticLat = 3.44398752157308;
  vtTarget.at(68).GalacticLat = 3.36218902541756;
  vtTarget.at(69).GalacticLat = 2.58088086166479;
  vtTarget.at(70).GalacticLat = -4.5225032741205;
  vtTarget.at(71).GalacticLat = -4.52502982145663;
  vtTarget.at(72).GalacticLat = -4.71802200241032;
  vtTarget.at(73).GalacticLat = -4.57523007937948;
  vtTarget.at(74).GalacticLat = -4.89263769586802;
  vtTarget.at(75).GalacticLat = -4.90269633150767;
  vtTarget.at(76).GalacticLat = -5.03890867670439;
  vtTarget.at(77).GalacticLat = -4.44498252241229;
  vtTarget.at(78).GalacticLat = -4.10271565732342;
  vtTarget.at(79).GalacticLat = -3.98535393107933;
  vtTarget.at(80).GalacticLat = -4.07910050208569;
  vtTarget.at(81).GalacticLat = -4.48699629387507;
  vtTarget.at(82).GalacticLat = -4.50367411651916;
  vtTarget.at(83).GalacticLat = -4.4958215461838;
  vtTarget.at(84).GalacticLat = -4.46899845081673;
  vtTarget.at(85).GalacticLat = -4.6961839270075;
  vtTarget.at(86).GalacticLat = -1.60598007878222;
  vtTarget.at(87).GalacticLat = -2.08528024930809;
  vtTarget.at(88).GalacticLat = -2.47121590996449;
  vtTarget.at(89).GalacticLat = -2.55353347013114;
  vtTarget.at(90).GalacticLat = -2.58849693181764;
  vtTarget.at(91).GalacticLat = -2.68968399884634;
  vtTarget.at(92).GalacticLat = -2.54296915848866;
  vtTarget.at(93).GalacticLat = -2.51622132813139;
  vtTarget.at(94).GalacticLat = -2.2329200001534;
  vtTarget.at(95).GalacticLat = -1.998445645641;
  vtTarget.at(96).GalacticLat = -2.87628155257424;
  vtTarget.at(97).GalacticLat = -2.49760643094393;
  vtTarget.at(98).GalacticLat = -2.602120870781;
  vtTarget.at(99).GalacticLat = -2.67129680154375;
  vtTarget.at(100).GalacticLat = -2.96335065294203;
  vtTarget.at(101).GalacticLat = -3.06358177381986;
  vtTarget.at(102).GalacticLat = -3.15634387669953;
  vtTarget.at(103).GalacticLat = -3.19687039511954;
  vtTarget.at(104).GalacticLat = -2.6596790630565;
  vtTarget.at(105).GalacticLat = -3.20687809529114;
  vtTarget.at(106).GalacticLat = -3.25675233612185;
  vtTarget.at(107).GalacticLat = -3.45869741813397;
  vtTarget.at(108).GalacticLat = -3.42751648283025;
  vtTarget.at(109).GalacticLat = -3.440373658107;
  vtTarget.at(110).GalacticLat = -3.45922996281465;
  vtTarget.at(111).GalacticLat = -1.05695242613923;
  vtTarget.at(112).GalacticLat = -1.45940488028191;
  vtTarget.at(113).GalacticLat = -1.15634688863114;
  vtTarget.at(114).GalacticLat = -1.52405088371641;
  vtTarget.at(115).GalacticLat = -1.34831554115572;
  vtTarget.at(116).GalacticLat = -1.46551009684022;
  vtTarget.at(117).GalacticLat = -1.53582452935616;
  vtTarget.at(118).GalacticLat = -1.611678895855;
  vtTarget.at(119).GalacticLat = -1.58081641149219;
  vtTarget.at(120).GalacticLat = -1.57212295019417;
  vtTarget.at(121).GalacticLat = -1.57977908351471;
  vtTarget.at(122).GalacticLat = -1.69939591859002;
  vtTarget.at(123).GalacticLat = -2.10001116405206;
  vtTarget.at(124).GalacticLat = -2.24347416350754;
  vtTarget.at(125).GalacticLat = -1.5960077189141;
  vtTarget.at(126).GalacticLat = -1.59718942451692;
  vtTarget.at(127).GalacticLat = -1.81076732489168;
  vtTarget.at(128).GalacticLat = -1.82517627976682;
  vtTarget.at(129).GalacticLat = -2.03427560593076;
  vtTarget.at(130).GalacticLat = -2.18237015614993;
  vtTarget.at(131).GalacticLat = -2.18929341927279;
  vtTarget.at(132).GalacticLat = -2.3174695835245;
  vtTarget.at(133).GalacticLat = -2.4518535185385;
  vtTarget.at(134).GalacticLat = -2.31431941415948;
  vtTarget.at(135).GalacticLat = -2.50474344502796;
  vtTarget.at(136).GalacticLat = -1.85944733831122;
  vtTarget.at(137).GalacticLat = -1.89326171979058;
  vtTarget.at(138).GalacticLat = -2.07752394292432;
  vtTarget.at(139).GalacticLat = -2.16093971289234;
  vtTarget.at(140).GalacticLat = -2.20156959031359;
  vtTarget.at(141).GalacticLat = -2.36224630103857;
  vtTarget.at(142).GalacticLat = -2.51546157226366;
  vtTarget.at(143).GalacticLat = -2.67413409689345;
  vtTarget.at(144).GalacticLat = -2.70887621635884;
  vtTarget.at(145).GalacticLat = -2.55963096750932;
  vtTarget.at(146).GalacticLat = -2.84339586392378;
  vtTarget.at(147).GalacticLat = -2.93644950077995;
  vtTarget.at(148).GalacticLat = -2.8400536186164;
  vtTarget.at(149).GalacticLat = -2.70312428876841;
  vtTarget.at(150).GalacticLat = -2.67373006581918;
  vtTarget.at(151).GalacticLat = -2.77536722872971;
  vtTarget.at(152).GalacticLat = -2.86674078822535;
  vtTarget.at(153).GalacticLat = -2.98404003937461;
  vtTarget.at(154).GalacticLat = -2.87973787684266;
  vtTarget.at(155).GalacticLat = -3.08740137601082;
  vtTarget.at(156).GalacticLat = -4.02100840768582;
  vtTarget.at(157).GalacticLat = -3.97462547114804;
  vtTarget.at(158).GalacticLat = -4.09062232540539;
  vtTarget.at(159).GalacticLat = -4.29012045790216;
  vtTarget.at(160).GalacticLat = -4.88210061054335;
  vtTarget.at(161).GalacticLat = -3.42207121493332;
  vtTarget.at(162).GalacticLat = -3.33867249972033;
  vtTarget.at(163).GalacticLat = -3.30258314337717;
  vtTarget.at(164).GalacticLat = -3.33589047779005;
  vtTarget.at(165).GalacticLat = -3.50637195219708;
  vtTarget.at(166).GalacticLat = -3.50045251634029;
  vtTarget.at(167).GalacticLat = -3.34127925750796;
  vtTarget.at(168).GalacticLat = -3.52093202301881;
  vtTarget.at(169).GalacticLat = -3.53337193731272;
  vtTarget.at(170).GalacticLat = -3.57279721070708;
  vtTarget.at(171).GalacticLat = -3.75282192780548;
  vtTarget.at(172).GalacticLat = -3.50412696295568;
  vtTarget.at(173).GalacticLat = -3.53159515506992;
  vtTarget.at(174).GalacticLat = -3.53554257802979;
  vtTarget.at(175).GalacticLat = -3.744255281255;
  vtTarget.at(176).GalacticLat = -3.89027737309466;
  vtTarget.at(177).GalacticLat = -3.54099704569742;
  vtTarget.at(178).GalacticLat = -4.18372706672996;
  vtTarget.at(179).GalacticLat = -4.80143870747463;
  vtTarget.at(180).GalacticLat = -4.57800274929379;
  vtTarget.at(181).GalacticLat = -2.84078575428615;
  vtTarget.at(182).GalacticLat = -2.88189125547237;
  vtTarget.at(183).GalacticLat = -3.21088725878944;
  vtTarget.at(184).GalacticLat = -4.02806462654408;
  vtTarget.at(185).GalacticLat = -4.3903607306228;
  vtTarget.at(186).GalacticLat = -3.01345957392618;
  vtTarget.at(187).GalacticLat = -3.13758390796338;
  vtTarget.at(188).GalacticLat = -3.00219364080696;
  vtTarget.at(189).GalacticLat = -4.15044629066694;
  vtTarget.at(190).GalacticLat = 1.26963536227204;
  vtTarget.at(191).GalacticLat = 1.41430552486818;
  vtTarget.at(192).GalacticLat = 1.56336067582695;
  vtTarget.at(193).GalacticLat = 1.16213502795027;
  vtTarget.at(194).GalacticLat = 2.06931716913102;
  vtTarget.at(195).GalacticLat = 2.01961290406855;
  vtTarget.at(196).GalacticLat = 1.70802885741178;
  vtTarget.at(197).GalacticLat = 1.42407376131587;
  vtTarget.at(198).GalacticLat = 2.91781324631301;
  vtTarget.at(199).GalacticLat = 2.83095865054314;
  vtTarget.at(200).GalacticLat = 4.37163833054845;
  vtTarget.at(201).GalacticLat = 4.33109474248165;
  vtTarget.at(202).GalacticLat = 3.20730192030543;
  vtTarget.at(203).GalacticLat = -4.76343148265208;
  vtTarget.at(204).GalacticLat = -4.29612461915472;
  vtTarget.at(205).GalacticLat = -4.47129067903013;
  vtTarget.at(206).GalacticLat = -0.764689566906341;
  vtTarget.at(207).GalacticLat = -2.11839680031662;
  vtTarget.at(208).GalacticLat = 1.33800367443259;
  vtTarget.at(209).GalacticLat = 3.5156525925707;
  vtTarget.at(210).GalacticLat = 4.91614261742739;
  vtTarget.at(211).GalacticLat = 3.36104347132441;
  vtTarget.at(212).GalacticLat = 3.21612777777539;
  vtTarget.at(213).GalacticLat = 2.72358695517585;
  vtTarget.at(214).GalacticLat = 3.7775367899998;
  vtTarget.at(215).GalacticLat = -1.30854667019713;
  vtTarget.at(216).GalacticLat = -2.39458163643293;
  vtTarget.at(217).GalacticLat = -3.7124175159804;
  vtTarget.at(218).GalacticLat = -0.871371843838351;
  vtTarget.at(219).GalacticLat = -0.985174256677115;
  vtTarget.at(220).GalacticLat = -1.1136598319426;
  vtTarget.at(221).GalacticLat = -1.44216107941819;
  vtTarget.at(222).GalacticLat = -3.3409789017383;
  vtTarget.at(223).GalacticLat = -3.55397370434572;
  vtTarget.at(224).GalacticLat = -3.66901635329464;
  vtTarget.at(225).GalacticLat = -1.10558388162007;
  vtTarget.at(226).GalacticLat = -1.46917403261551;
  vtTarget.at(227).GalacticLat = -3.63521970785518;
  vtTarget.at(228).GalacticLat = -4.10034037525346;
  vtTarget.at(229).GalacticLat = -4.29183221932043;
  vtTarget.at(230).GalacticLat = -1.03905965735693;
  vtTarget.at(231).GalacticLat = -1.03564683819331;
  vtTarget.at(232).GalacticLat = -1.43665362699968;
  vtTarget.at(233).GalacticLat = -1.55223447931697;
  vtTarget.at(234).GalacticLat = -1.54292976317559;
  vtTarget.at(235).GalacticLat = -1.71207311633674;
  vtTarget.at(236).GalacticLat = -3.34765365009339;
  vtTarget.at(237).GalacticLat = -1.38865883937359;
  vtTarget.at(238).GalacticLat = -1.81480039654679;
  vtTarget.at(239).GalacticLat = -1.8125762127201;
  vtTarget.at(240).GalacticLat = -1.85274527260641;
  vtTarget.at(241).GalacticLat = -1.84884765834647;
  vtTarget.at(242).GalacticLat = -2.31758931567376;
  vtTarget.at(243).GalacticLat = -1.62730338729531;
  vtTarget.at(244).GalacticLat = -2.00482066280152;
  vtTarget.at(245).GalacticLat = -2.05677397358753;
  vtTarget.at(246).GalacticLat = -2.06739026408889;
  vtTarget.at(247).GalacticLat = -2.23426723160726;
  vtTarget.at(248).GalacticLat = -2.42126380973104;
  vtTarget.at(249).GalacticLat = -2.46785006910472;
  vtTarget.at(250).GalacticLat = -3.221888122676;
  vtTarget.at(251).GalacticLat = -3.24256153773936;
  vtTarget.at(252).GalacticLat = -2.45144854719141;
  vtTarget.at(253).GalacticLat = -2.65538698915252;
  vtTarget.at(254).GalacticLat = -2.65824474391319;
  vtTarget.at(255).GalacticLat = -2.6011825751852;
  vtTarget.at(256).GalacticLat = -2.66385111197698;
  vtTarget.at(257).GalacticLat = -2.67467928806562;
  vtTarget.at(258).GalacticLat = -2.89942184990679;
  vtTarget.at(259).GalacticLat = -2.84780045221562;
  vtTarget.at(260).GalacticLat = -2.78290432124273;
  vtTarget.at(261).GalacticLat = -4.51451645893655;
  vtTarget.at(262).GalacticLat = -4.51693308766455;
  vtTarget.at(263).GalacticLat = -3.03704265855366;
  vtTarget.at(264).GalacticLat = -3.2672706444869;
  vtTarget.at(265).GalacticLat = -3.28090675569896;
  vtTarget.at(266).GalacticLat = -2.41352404704863;
  vtTarget.at(267).GalacticLat = -2.46951690018393;
  vtTarget.at(268).GalacticLat = -3.16741530274695;
  vtTarget.at(269).GalacticLat = -3.22090721262077;
  vtTarget.at(270).GalacticLat = -3.28836301392095;
  vtTarget.at(271).GalacticLat = -3.26369954661861;
  vtTarget.at(272).GalacticLat = -2.64935585536232;
  vtTarget.at(273).GalacticLat = -3.20663860511452;
  vtTarget.at(274).GalacticLat = -3.19934376694093;
  vtTarget.at(275).GalacticLat = -3.62684802624993;
  vtTarget.at(276).GalacticLat = -4.78855288386557;
  vtTarget.at(277).GalacticLat = -0.228631451529957;
  vtTarget.at(278).GalacticLat = 4.02734561873738;
  vtTarget.at(279).GalacticLat = 1.805847046558;
  vtTarget.at(280).GalacticLat = -3.21516787944223;
  vtTarget.at(281).GalacticLat = -2.33818903963567;
  vtTarget.at(282).GalacticLat = -4.62401132858022;
  vtTarget.at(283).GalacticLat = -3.79766833912534;
  vtTarget.at(284).GalacticLat = -0.783459926484891;
  vtTarget.at(285).GalacticLat = -2.03718345040172;
  vtTarget.at(286).GalacticLat = 2.86087876904375;
  vtTarget.at(287).GalacticLat = -4.41667322843999;
  vtTarget.at(288).GalacticLat = 1.76462746861067;
  vtTarget.at(289).GalacticLat = -3.09350831447548;
  vtTarget.at(290).GalacticLat = 3.96272774925538;
  vtTarget.at(291).GalacticLat = -4.75165500405651;
  vtTarget.at(292).GalacticLat = -3.72562809518669;
  vtTarget.at(293).GalacticLat = -2.91072752473508;
  vtTarget.at(294).GalacticLat = 1.21286444788372;
  vtTarget.at(295).GalacticLat = 1.86140846949933;
  vtTarget.at(296).GalacticLat = -1.0674617123102;
  vtTarget.at(297).GalacticLat = -1.21944581757473;
  vtTarget.at(298).GalacticLat = -0.509948165778625;
  vtTarget.at(299).GalacticLat = -4.62011243235467;
  vtTarget.at(300).GalacticLat = -3.27180477688659;
  vtTarget.at(301).GalacticLat = -0.0250097651060424;
  vtTarget.at(302).GalacticLat = -0.0461210248267241;
  vtTarget.at(303).GalacticLat = 2.02508873348652;
  vtTarget.at(304).GalacticLat = -4.66679559726563;
  vtTarget.at(305).GalacticLat = 0.332322552114582;
  vtTarget.at(306).GalacticLat = -1.3835151700082;
  vtTarget.at(307).GalacticLat = -0.0882055888885365;
  vtTarget.at(308).GalacticLat = 2.7014605960244;
  vtTarget.at(309).GalacticLat = -3.96188837710361;
  vtTarget.at(310).GalacticLat = -3.62686698966618;
  vtTarget.at(311).GalacticLat = 0.853067853030397;
  vtTarget.at(312).GalacticLat = -1.87218084970269;
  vtTarget.at(313).GalacticLat = -1.39259216194154;
  vtTarget.at(314).GalacticLat = 1.32084062011572;
  vtTarget.at(315).GalacticLat = -4.58295955607195;
  vtTarget.at(316).GalacticLat = -3.71490320038522;
  vtTarget.at(317).GalacticLat = -2.71470379274134;
  vtTarget.at(318).GalacticLat = -3.87953308771147;
  vtTarget.at(319).GalacticLat = -2.07063018553748;
  vtTarget.at(320).GalacticLat = 4.8883661586959;
  vtTarget.at(321).GalacticLat = 3.91350061239352;
  vtTarget.at(322).GalacticLat = -1.27952400107486;
  vtTarget.at(323).GalacticLat = -4.42372638378314;
  vtTarget.at(324).GalacticLat = -1.4539307193518;
  vtTarget.at(325).GalacticLat = 3.96212004026235;
  vtTarget.at(326).GalacticLat = 1.20402780290139;
  vtTarget.at(327).GalacticLat = 3.13466391011242;
  vtTarget.at(328).GalacticLat = 1.31999147301951;
  vtTarget.at(329).GalacticLat = 4.28703410145045;
  vtTarget.at(330).GalacticLat = 4.08493660369583;
  vtTarget.at(331).GalacticLat = 4.92090615784888;
  vtTarget.at(332).GalacticLat = 3.38959423350969;
  vtTarget.at(333).GalacticLat = 4.29559544077394;
  vtTarget.at(334).GalacticLat = 3.54177016828128;
  vtTarget.at(335).GalacticLat = -3.66724895091302;
  vtTarget.at(336).GalacticLat = 4.88790461286628;
  vtTarget.at(337).GalacticLat = 1.25503771379294;
  vtTarget.at(338).GalacticLat = -2.58475883456784;
  vtTarget.at(339).GalacticLat = -4.21380293161863;
  vtTarget.at(340).GalacticLat = 0.239258381357445;
  vtTarget.at(341).GalacticLat = 1.43423183929664;
  vtTarget.at(342).GalacticLat = 2.82061667802415;
  vtTarget.at(343).GalacticLat = 4.46246309705355;
  vtTarget.at(344).GalacticLat = -1.83861002140009;
  vtTarget.at(345).GalacticLat = 4.54333856305584;
  vtTarget.at(346).GalacticLat = 2.93807949339195;
  vtTarget.at(347).GalacticLat = -3.97187046316319;
  vtTarget.at(348).GalacticLat = -0.423912437792426;
  vtTarget.at(349).GalacticLat = 3.76632862317367;
  vtTarget.at(350).GalacticLat = -1.50130270277275;
  vtTarget.at(351).GalacticLat = -3.12369620690709;
  vtTarget.at(352).GalacticLat = -2.00575751900137;
  vtTarget.at(353).GalacticLat = -3.11929916872297;
  vtTarget.at(354).GalacticLat = 1.81116469659451;
  vtTarget.at(355).GalacticLat = 2.54843392402452;
  vtTarget.at(356).GalacticLat = -4.76206614044804;
  vtTarget.at(357).GalacticLat = -2.49590991853435;
  vtTarget.at(358).GalacticLat = 3.31557763901895;
  vtTarget.at(359).GalacticLat = -1.55110874238914;
  vtTarget.at(360).GalacticLat = -1.86871707283748;
  vtTarget.at(361).GalacticLat = -0.428832634690652;
  vtTarget.at(362).GalacticLat = 1.02626549232564;
  vtTarget.at(363).GalacticLat = 4.32590909424328;
  vtTarget.at(364).GalacticLat = -4.12115230925706;
  vtTarget.at(365).GalacticLat = 0.432214547134278;
  vtTarget.at(366).GalacticLat = -3.08145667988468;
  vtTarget.at(367).GalacticLat = 1.06348245683737;
  vtTarget.at(368).GalacticLat = 3.57379087272255;
  vtTarget.at(369).GalacticLat = -1.84729196789263;
  vtTarget.at(370).GalacticLat = 2.45873210598963;
  vtTarget.at(371).GalacticLat = -4.13857145349229;
  vtTarget.at(372).GalacticLat = -3.04044027202318;
  vtTarget.at(373).GalacticLat = -4.74050866607999;
  vtTarget.at(374).GalacticLat = -3.96049903958301;
  vtTarget.at(375).GalacticLat = -4.29922934074024;
  vtTarget.at(376).GalacticLat = 0.769097039838832;
  vtTarget.at(377).GalacticLat = -1.75403899163259;
  vtTarget.at(378).GalacticLat = 2.39039734735146;
  vtTarget.at(379).GalacticLat = -2.95781280685511;
  vtTarget.at(380).GalacticLat = 4.4145499183629;
  vtTarget.at(381).GalacticLat = -2.41305854652083;
  vtTarget.at(382).GalacticLat = 4.65215234036928;
  vtTarget.at(383).GalacticLat = -3.31887703057352;
  vtTarget.at(384).GalacticLat = -3.89623864176841;
  vtTarget.at(385).GalacticLat = -0.920398290939254;
  vtTarget.at(386).GalacticLat = -4.37130327566451;
  vtTarget.at(387).GalacticLat = 3.77540984132983;
  vtTarget.at(388).GalacticLat = 0.290423318838259;
  vtTarget.at(389).GalacticLat = -3.38406569150685;
  vtTarget.at(390).GalacticLat = -4.52743334030108;
  vtTarget.at(391).GalacticLat = -0.175139474104711;
  vtTarget.at(392).GalacticLat = -1.94437975109958;
  vtTarget.at(393).GalacticLat = 2.24492534197946;
  vtTarget.at(394).GalacticLat = 2.56693756120307;
  vtTarget.at(395).GalacticLat = -0.69959087617779;
  vtTarget.at(396).GalacticLat = -4.99367985715302;
  vtTarget.at(397).GalacticLat = 4.87913384547885;
  vtTarget.at(398).GalacticLat = -0.578061738309614;
  vtTarget.at(399).GalacticLat = 3.46312190471298;
  vtTarget.at(400).GalacticLat = -0.724431058060214;
  vtTarget.at(401).GalacticLat = -0.678053955706647;
  vtTarget.at(402).GalacticLat = 1.68529556006996;
  vtTarget.at(403).GalacticLat = -1.5409814636188;
  vtTarget.at(404).GalacticLat = 0.111902043095212;
  vtTarget.at(405).GalacticLat = 1.48484081339841;
  vtTarget.at(406).GalacticLat = 3.21439332535815;
  vtTarget.at(407).GalacticLat = -1.29130056469273;
  vtTarget.at(408).GalacticLat = 0.634977576456009;
  vtTarget.at(409).GalacticLat = -4.3608916379065;
  vtTarget.at(410).GalacticLat = -2.53121658375012;
  vtTarget.at(411).GalacticLat = -1.34170576739395;
  vtTarget.at(412).GalacticLat = 1.65484024051216;
  vtTarget.at(413).GalacticLat = -4.23630090057162;
  vtTarget.at(414).GalacticLat = 0.321597027787369;
  vtTarget.at(415).GalacticLat = 0.308779973433975;
  vtTarget.at(416).GalacticLat = 3.74814269416808;
  vtTarget.at(417).GalacticLat = 3.39008610164271;
  vtTarget.at(418).GalacticLat = 2.58648495918829;
  vtTarget.at(419).GalacticLat = -2.2258774317467;
  vtTarget.at(420).GalacticLat = -1.83043913560324;
  vtTarget.at(421).GalacticLat = -3.61061447212958;
  vtTarget.at(422).GalacticLat = 0.21137918681077;
  vtTarget.at(423).GalacticLat = 1.13691477510892;
  vtTarget.at(424).GalacticLat = 4.48141984474941;
  vtTarget.at(425).GalacticLat = 3.40668477909254;
  vtTarget.at(426).GalacticLat = -0.657726233148388;
  vtTarget.at(427).GalacticLat = 3.22365724235316;
  vtTarget.at(428).GalacticLat = -4.83774059232013;
  vtTarget.at(429).GalacticLat = -1.17558893704471;
  vtTarget.at(430).GalacticLat = -4.80012111600332;
  vtTarget.at(431).GalacticLat = 1.56348157235274;
  vtTarget.at(432).GalacticLat = 3.2522433332774;
  vtTarget.at(433).GalacticLat = -0.0636927361735147;
  vtTarget.at(434).GalacticLat = -2.94672382136354;
  vtTarget.at(435).GalacticLat = -4.09079474588811;
  vtTarget.at(436).GalacticLat = 4.58529697740816;
  vtTarget.at(437).GalacticLat = -0.126916387883092;
  vtTarget.at(438).GalacticLat = 2.0199577696448;
  vtTarget.at(439).GalacticLat = 0.0188532087545479;
  vtTarget.at(440).GalacticLat = 0.82324475160102;
  vtTarget.at(441).GalacticLat = 3.32807291678724;
  vtTarget.at(442).GalacticLat = -3.9266662413024;
  vtTarget.at(443).GalacticLat = 2.6394562670268;
  vtTarget.at(444).GalacticLat = -0.507126892594889;
  vtTarget.at(445).GalacticLat = -2.90025885962825;
  vtTarget.at(446).GalacticLat = 0.10001190358653;
  vtTarget.at(447).GalacticLat = 0.648784778770331;
  vtTarget.at(448).GalacticLat = 3.59184846834268;
  vtTarget.at(449).GalacticLat = 2.13363037285712;
  vtTarget.at(450).GalacticLat = 4.58329272453103;
  vtTarget.at(451).GalacticLat = 1.13366254674649;
  vtTarget.at(452).GalacticLat = 4.88564827612275;
  vtTarget.at(453).GalacticLat = 3.74694652463736;
  vtTarget.at(454).GalacticLat = 4.39147862161348;
  vtTarget.at(455).GalacticLat = 1.83430241326866;
  vtTarget.at(456).GalacticLat = -4.76587771143152;
  vtTarget.at(457).GalacticLat = -3.15482547825052;
  vtTarget.at(458).GalacticLat = 1.88344300130143;
  vtTarget.at(459).GalacticLat = 1.04377564684204;
  vtTarget.at(460).GalacticLat = -1.2338391559786;
  vtTarget.at(461).GalacticLat = 2.12543180411929;
  vtTarget.at(462).GalacticLat = -3.63840297218842;
  vtTarget.at(463).GalacticLat = -0.277551305376229;
  vtTarget.at(464).GalacticLat = 1.34290303468472;
  vtTarget.at(465).GalacticLat = -1.7285712458662;
  vtTarget.at(466).GalacticLat = 0.339872582760612;
  vtTarget.at(467).GalacticLat = 1.83302169469654;
  vtTarget.at(468).GalacticLat = 4.11455178446355;
  vtTarget.at(469).GalacticLat = 4.72572767339935;
  vtTarget.at(470).GalacticLat = -3.31632120336536;
  vtTarget.at(471).GalacticLat = -3.38054159929401;
  vtTarget.at(472).GalacticLat = -0.425041976939943;
  vtTarget.at(473).GalacticLat = -0.0881621348072539;
  vtTarget.at(474).GalacticLat = 4.57352808012363;
  vtTarget.at(475).GalacticLat = -2.5086002473461;
  vtTarget.at(476).GalacticLat = 4.23349117187282;
  vtTarget.at(477).GalacticLat = 4.42264730201683;
  vtTarget.at(478).GalacticLat = 1.64009537506927;
  vtTarget.at(479).GalacticLat = 2.88557198643927;
  vtTarget.at(480).GalacticLat = -2.0355127257548;
  vtTarget.at(481).GalacticLat = -0.188639382977767;
  vtTarget.at(482).GalacticLat = -1.93022432979582;
  vtTarget.at(483).GalacticLat = -4.32627436160638;
  vtTarget.at(484).GalacticLat = 3.3954619564112;
  vtTarget.at(485).GalacticLat = -2.38493369581981;
  vtTarget.at(486).GalacticLat = -2.45507536373211;
  vtTarget.at(487).GalacticLat = 3.96441731231128;
  vtTarget.at(488).GalacticLat = 2.7699880818681;
  vtTarget.at(489).GalacticLat = 2.33065769674582;
  vtTarget.at(490).GalacticLat = -3.48636329690553;
  vtTarget.at(491).GalacticLat = -2.21619792881211;
  vtTarget.at(492).GalacticLat = 1.50459345802662;
  vtTarget.at(493).GalacticLat = 2.73148639546358;
  vtTarget.at(494).GalacticLat = -2.54123547419735;
  vtTarget.at(495).GalacticLat = 2.67341384256252;
  vtTarget.at(496).GalacticLat = 4.15842035307874;
  vtTarget.at(497).GalacticLat = -1.23529224437823;
  vtTarget.at(498).GalacticLat = -3.67492058543133;
  vtTarget.at(499).GalacticLat = -0.518604497489294;
  vtTarget.at(500).GalacticLat = 2.27042652295165;
  vtTarget.at(501).GalacticLat = 3.45443784110587;
  vtTarget.at(502).GalacticLat = 0.424232560732367;
  vtTarget.at(503).GalacticLat = 1.15782096989085;
  vtTarget.at(504).GalacticLat = 1.26250183610036;
  vtTarget.at(505).GalacticLat = -0.247977831685016;
  vtTarget.at(506).GalacticLat = 4.27103328874135;
  vtTarget.at(507).GalacticLat = 1.11273361550866;
  vtTarget.at(508).GalacticLat = -2.81190792175084;
  vtTarget.at(509).GalacticLat = -3.7642283428376;
  vtTarget.at(510).GalacticLat = -0.128074469321575;
  vtTarget.at(511).GalacticLat = -2.6154278504668;
  vtTarget.at(512).GalacticLat = -2.18535982658246;
  vtTarget.at(513).GalacticLat = -0.102391439211188;
  vtTarget.at(514).GalacticLat = -4.93302302811991;
  vtTarget.at(515).GalacticLat = -4.63697364191794;
  vtTarget.at(516).GalacticLat = 4.31675718403067;
  vtTarget.at(517).GalacticLat = 0.0389250530517007;
  vtTarget.at(518).GalacticLat = 3.91757225812212;
  vtTarget.at(519).GalacticLat = -2.57707450129068;
  vtTarget.at(520).GalacticLat = 4.80249451740616;
  vtTarget.at(521).GalacticLat = -3.36872774005767;
  vtTarget.at(522).GalacticLat = -1.55995677227624;
  vtTarget.at(523).GalacticLat = 1.03834816447245;
  vtTarget.at(524).GalacticLat = 3.23234647474367;
  vtTarget.at(525).GalacticLat = -3.92993193334511;
  vtTarget.at(526).GalacticLat = 4.10943808663849;
  vtTarget.at(527).GalacticLat = 4.99738096958076;
  vtTarget.at(528).GalacticLat = 4.49722843222794;
  vtTarget.at(529).GalacticLat = 0.227156937170614;
  vtTarget.at(530).GalacticLat = -4.40304072162234;
  vtTarget.at(531).GalacticLat = -0.280657027449463;
  vtTarget.at(532).GalacticLat = 4.44048600840138;
  vtTarget.at(533).GalacticLat = 1.60128055390698;
  vtTarget.at(534).GalacticLat = -2.16467957626828;
  vtTarget.at(535).GalacticLat = 2.84853853075362;
  vtTarget.at(536).GalacticLat = -1.83494317567094;
  vtTarget.at(537).GalacticLat = -3.93679121107654;
  vtTarget.at(538).GalacticLat = -4.12817782951729;
  vtTarget.at(539).GalacticLat = 3.84450695878445;
  vtTarget.at(540).GalacticLat = 3.35246038030585;
  vtTarget.at(541).GalacticLat = -0.575351988143405;
  vtTarget.at(542).GalacticLat = 0.458723292000011;
  vtTarget.at(543).GalacticLat = 2.47840094939029;
  vtTarget.at(544).GalacticLat = 1.65809016565331;
  vtTarget.at(545).GalacticLat = -0.875054935088705;
  vtTarget.at(546).GalacticLat = -2.81128063062686;
  vtTarget.at(547).GalacticLat = 3.28834902021459;
  vtTarget.at(548).GalacticLat = 4.98772635452227;
  vtTarget.at(549).GalacticLat = 2.02853549810098;
  vtTarget.at(550).GalacticLat = -0.282006256409345;
  vtTarget.at(551).GalacticLat = -0.0644674699339154;
  vtTarget.at(552).GalacticLat = -4.36088696951409;
  vtTarget.at(553).GalacticLat = -2.03639233632618;
  vtTarget.at(554).GalacticLat = 0.541853413991525;
  vtTarget.at(555).GalacticLat = -4.80441232089296;
  vtTarget.at(556).GalacticLat = -4.00000023248682;
  vtTarget.at(557).GalacticLat = 0.355137020696819;
  vtTarget.at(558).GalacticLat = -2.94119948416467;
  vtTarget.at(559).GalacticLat = 1.62478014653675;
  vtTarget.at(560).GalacticLat = 0.215184236781893;
  vtTarget.at(561).GalacticLat = -2.66838021025178;
  vtTarget.at(562).GalacticLat = -1.00599383248521;
  vtTarget.at(563).GalacticLat = 3.20346526749988;
  vtTarget.at(564).GalacticLat = -1.10683282987754;
  vtTarget.at(565).GalacticLat = 1.19370327024349;
  vtTarget.at(566).GalacticLat = 4.27237708110692;
  vtTarget.at(567).GalacticLat = 3.96303051103363;
  vtTarget.at(568).GalacticLat = -4.83864094300033;
  vtTarget.at(569).GalacticLat = -0.0301661735325487;
  vtTarget.at(570).GalacticLat = -0.289431117504583;
  vtTarget.at(571).GalacticLat = -0.113019070479876;
  vtTarget.at(572).GalacticLat = 0.523630397109995;
  vtTarget.at(573).GalacticLat = -1.50211497320706;
  vtTarget.at(574).GalacticLat = -3.73573383839527;
  vtTarget.at(575).GalacticLat = 0.550317186301309;
  vtTarget.at(576).GalacticLat = -0.581179238609458;
  vtTarget.at(577).GalacticLat = -4.99292288035404;
  vtTarget.at(578).GalacticLat = 1.13667353999718;
  vtTarget.at(579).GalacticLat = 3.120557874407;
  vtTarget.at(580).GalacticLat = -0.271998089665081;
  vtTarget.at(581).GalacticLat = -1.7802337542738;
  vtTarget.at(582).GalacticLat = -3.04616254622611;
  vtTarget.at(583).GalacticLat = -4.68770709644682;
  vtTarget.at(584).GalacticLat = 1.35022760195378;
  vtTarget.at(585).GalacticLat = -4.8844522568365;
  vtTarget.at(586).GalacticLat = 1.36158604353108;
  vtTarget.at(587).GalacticLat = -2.96767128552852;
  vtTarget.at(588).GalacticLat = -0.697211004831672;
  vtTarget.at(589).GalacticLat = 0.689501620817285;
  vtTarget.at(590).GalacticLat = 4.17710500038929;
  vtTarget.at(591).GalacticLat = 4.11778552213937;
  vtTarget.at(592).GalacticLat = 3.32742670196268;
  vtTarget.at(593).GalacticLat = -4.85080329189086;
  vtTarget.at(594).GalacticLat = -2.57223757819463;
  vtTarget.at(595).GalacticLat = -2.3816538625327;
  vtTarget.at(596).GalacticLat = -4.571334856243;
  vtTarget.at(597).GalacticLat = -0.449851390592617;
  vtTarget.at(598).GalacticLat = -2.28821870737145;
  vtTarget.at(599).GalacticLat = -1.26697190894393;
  vtTarget.at(600).GalacticLat = -2.97412560891041;
  vtTarget.at(601).GalacticLat = -3.42161823916292;
  vtTarget.at(602).GalacticLat = -3.50536842726825;
  vtTarget.at(603).GalacticLat = 1.7019337757872;
  vtTarget.at(604).GalacticLat = 4.47419050057004;
  vtTarget.at(605).GalacticLat = 2.33024763262161;
  vtTarget.at(606).GalacticLat = 2.89317065456875;
  vtTarget.at(607).GalacticLat = -0.325589064961189;
  vtTarget.at(608).GalacticLat = -1.84012186881659;
  vtTarget.at(609).GalacticLat = -2.86945856898615;
  vtTarget.at(610).GalacticLat = 3.50533960453678;
  vtTarget.at(611).GalacticLat = -3.04668998532377;
  vtTarget.at(612).GalacticLat = -3.94865870713174;
  vtTarget.at(613).GalacticLat = 1.70435569045757;
  vtTarget.at(614).GalacticLat = 4.90752657546489;
  vtTarget.at(615).GalacticLat = -3.23970724400438;
  vtTarget.at(616).GalacticLat = -2.05283953413486;
  vtTarget.at(617).GalacticLat = 2.31334957023106;
  vtTarget.at(618).GalacticLat = -0.186434034452732;
  vtTarget.at(619).GalacticLat = 3.59186399274413;
  vtTarget.at(620).GalacticLat = 3.49061024612865;
  vtTarget.at(621).GalacticLat = -3.09017904847732;
  vtTarget.at(622).GalacticLat = -0.222786995424901;
  vtTarget.at(623).GalacticLat = 4.56024988134992;
  vtTarget.at(624).GalacticLat = -0.524518804110524;
  vtTarget.at(625).GalacticLat = -2.15111035600393;
  vtTarget.at(626).GalacticLat = 4.31165734705275;
  vtTarget.at(627).GalacticLat = 3.5187770574125;
  vtTarget.at(628).GalacticLat = 4.1996130780544;
  vtTarget.at(629).GalacticLat = -4.82293974883047;
  vtTarget.at(630).GalacticLat = -4.20135194722752;
  vtTarget.at(631).GalacticLat = -4.08134731517227;
  vtTarget.at(632).GalacticLat = -4.94390408975988;
  vtTarget.at(633).GalacticLat = -1.72324187866354;
  vtTarget.at(634).GalacticLat = 4.19664910815862;
  vtTarget.at(635).GalacticLat = -3.14269212461179;
  vtTarget.at(636).GalacticLat = -2.2409899057876;
  vtTarget.at(637).GalacticLat = -4.38408565835956;
  vtTarget.at(638).GalacticLat = 1.9015224170538;
  vtTarget.at(639).GalacticLat = -0.561693163074784;
  vtTarget.at(640).GalacticLat = 3.81392350819364;
  vtTarget.at(641).GalacticLat = 2.15327564196839;
  vtTarget.at(642).GalacticLat = -2.19822988111376;
  vtTarget.at(643).GalacticLat = -2.9356438033003;
  vtTarget.at(644).GalacticLat = 1.22360878115694;
  vtTarget.at(645).GalacticLat = -4.03247879782025;
  vtTarget.at(646).GalacticLat = -2.55724588760725;
  vtTarget.at(647).GalacticLat = 2.88370679223416;
  vtTarget.at(648).GalacticLat = -2.93640226257896;
  vtTarget.at(649).GalacticLat = 0.441651601523711;
  vtTarget.at(650).GalacticLat = -4.37416804610543;
  vtTarget.at(651).GalacticLat = -0.914088973289706;
  vtTarget.at(652).GalacticLat = -0.230717369811097;
  vtTarget.at(653).GalacticLat = -4.56094202962348;
  vtTarget.at(654).GalacticLat = -0.353941722722139;
  vtTarget.at(655).GalacticLat = 2.99912461554415;
  vtTarget.at(656).GalacticLat = 4.54537559502587;
  vtTarget.at(657).GalacticLat = 4.50700335898505;
  vtTarget.at(658).GalacticLat = -1.17215830490726;
  vtTarget.at(659).GalacticLat = 2.6527727661057;
  vtTarget.at(660).GalacticLat = 3.81346561235536;
  vtTarget.at(661).GalacticLat = 3.33756205883049;
  vtTarget.at(662).GalacticLat = -4.41903092486828;
  vtTarget.at(663).GalacticLat = 4.89834194281235;
  vtTarget.at(664).GalacticLat = -2.67131893021298;
  vtTarget.at(665).GalacticLat = 2.44697693512089;
  vtTarget.at(666).GalacticLat = -4.45470456610837;
  vtTarget.at(667).GalacticLat = 4.76358032651952;
  vtTarget.at(668).GalacticLat = -4.81587870072648;
  vtTarget.at(669).GalacticLat = 2.51746623250989;
  vtTarget.at(670).GalacticLat = 1.10133602723093;
  vtTarget.at(671).GalacticLat = -0.852291836879273;
  vtTarget.at(672).GalacticLat = -0.348390745824153;
  vtTarget.at(673).GalacticLat = -2.42781211628585;
  vtTarget.at(674).GalacticLat = -3.62287843344046;
  vtTarget.at(675).GalacticLat = 2.73084601085361;
  vtTarget.at(676).GalacticLat = -1.76177684075943;
  vtTarget.at(677).GalacticLat = -1.90195554062248;
  vtTarget.at(678).GalacticLat = 4.56685317421099;
  vtTarget.at(679).GalacticLat = 0.253603272668999;
  vtTarget.at(680).GalacticLat = -2.90644401842168;
  vtTarget.at(681).GalacticLat = 1.72802277471864;
  vtTarget.at(682).GalacticLat = 0.746819249933331;
  vtTarget.at(683).GalacticLat = -4.81396934510449;
  vtTarget.at(684).GalacticLat = -4.09617338804793;
  vtTarget.at(685).GalacticLat = -4.33055843657061;
  vtTarget.at(686).GalacticLat = 3.15934648191975;
  vtTarget.at(687).GalacticLat = 4.1463001468846;
  vtTarget.at(688).GalacticLat = 4.53973364967527;
  vtTarget.at(689).GalacticLat = 2.26460195155973;
  vtTarget.at(690).GalacticLat = -4.80882098659979;
  vtTarget.at(691).GalacticLat = -4.77909839119874;
  vtTarget.at(692).GalacticLat = 3.78376814396976;
  vtTarget.at(693).GalacticLat = -2.5362589141409;
  vtTarget.at(694).GalacticLat = -1.80665961477052;
  vtTarget.at(695).GalacticLat = 3.07181133484114;
  vtTarget.at(696).GalacticLat = 0.659789743573223;
  vtTarget.at(697).GalacticLat = -0.557061721089217;
  vtTarget.at(698).GalacticLat = -0.731280008052217;
  vtTarget.at(699).GalacticLat = -2.81892899064538;
  vtTarget.at(700).GalacticLat = 4.43917523189899;
  vtTarget.at(701).GalacticLat = -1.69205526749555;
  vtTarget.at(702).GalacticLat = 3.53297289144418;
  vtTarget.at(703).GalacticLat = 2.78691132121866;
  vtTarget.at(704).GalacticLat = 1.25532561984341;
  vtTarget.at(705).GalacticLat = 4.96952979590712;
  vtTarget.at(706).GalacticLat = 1.76397637780841;
  vtTarget.at(707).GalacticLat = -2.05646236790249;
  vtTarget.at(708).GalacticLat = 3.46555122978031;
  vtTarget.at(709).GalacticLat = 1.42155302631771;
  vtTarget.at(710).GalacticLat = -4.62897888784795;
  vtTarget.at(711).GalacticLat = 4.49856279241902;
  vtTarget.at(712).GalacticLat = -4.14793283735117;
  vtTarget.at(713).GalacticLat = -2.64441900520066;
  vtTarget.at(714).GalacticLat = -0.794004031006851;
  vtTarget.at(715).GalacticLat = 4.36541914557319;
  vtTarget.at(716).GalacticLat = -2.00787034445379;
  vtTarget.at(717).GalacticLat = -0.392548737589097;
  vtTarget.at(718).GalacticLat = 2.23076094651892;
  vtTarget.at(719).GalacticLat = 3.61851380502772;
  vtTarget.at(720).GalacticLat = -4.30880916770188;
  vtTarget.at(721).GalacticLat = 1.12252108594223;
  vtTarget.at(722).GalacticLat = -2.58070216142557;
  vtTarget.at(723).GalacticLat = 1.21127901690931;
  vtTarget.at(724).GalacticLat = -0.263706899608809;
  vtTarget.at(725).GalacticLat = -1.89288358643474;
  vtTarget.at(726).GalacticLat = -3.26153142298714;
  vtTarget.at(727).GalacticLat = -4.63664496045723;
  vtTarget.at(728).GalacticLat = 2.00434446174644;
  vtTarget.at(729).GalacticLat = -4.07618390105931;
  vtTarget.at(730).GalacticLat = -3.55174118154347;
  vtTarget.at(731).GalacticLat = 0.758085684719243;
  vtTarget.at(732).GalacticLat = 4.78886303242243;
  vtTarget.at(733).GalacticLat = 3.49641943043594;
  vtTarget.at(734).GalacticLat = 1.2876926494675;
  vtTarget.at(735).GalacticLat = 3.52268206473019;
  vtTarget.at(736).GalacticLat = 2.77964591824012;
  vtTarget.at(737).GalacticLat = 3.1807890107965;
  vtTarget.at(738).GalacticLat = 2.71719872742379;
  vtTarget.at(739).GalacticLat = -2.02089480983615;
  vtTarget.at(740).GalacticLat = -1.67315749807153;
  vtTarget.at(741).GalacticLat = -1.07522255761188;
  vtTarget.at(742).GalacticLat = 0.164957529994291;
  vtTarget.at(743).GalacticLat = 3.66495850920367;
  vtTarget.at(744).GalacticLat = 2.68191526506896;
  vtTarget.at(745).GalacticLat = -3.68901461185843;
  vtTarget.at(746).GalacticLat = 3.84687855857336;
  vtTarget.at(747).GalacticLat = 1.94596908654448;
  vtTarget.at(748).GalacticLat = -2.58119572730905;
  vtTarget.at(749).GalacticLat = 2.9482186431477;
  vtTarget.at(750).GalacticLat = -4.32393178324957;
  vtTarget.at(751).GalacticLat = 4.47759900232862;
  vtTarget.at(752).GalacticLat = -2.08833353955746;
  vtTarget.at(753).GalacticLat = -0.733758489117331;
  vtTarget.at(754).GalacticLat = -3.74345262103477;
  vtTarget.at(755).GalacticLat = -0.718916118912157;
  vtTarget.at(756).GalacticLat = 4.12724303008593;
  vtTarget.at(757).GalacticLat = -0.0962468807362988;
  vtTarget.at(758).GalacticLat = 1.95047389289411;
  vtTarget.at(759).GalacticLat = 1.67310271279509;
  vtTarget.at(760).GalacticLat = -0.292968813878545;
  vtTarget.at(761).GalacticLat = 2.96616168204266;
  vtTarget.at(762).GalacticLat = -0.523329170432635;
  vtTarget.at(763).GalacticLat = 2.22051975264733;
  vtTarget.at(764).GalacticLat = -2.22958733883392;
  vtTarget.at(765).GalacticLat = 3.9836687521383;
  vtTarget.at(766).GalacticLat = 3.91222664957901;
  vtTarget.at(767).GalacticLat = -1.99461191855888;
  vtTarget.at(768).GalacticLat = -2.2070297113153;
  vtTarget.at(769).GalacticLat = 3.48181203688036;
  vtTarget.at(770).GalacticLat = -2.44519091492765;
  vtTarget.at(771).GalacticLat = -1.3863909642622;
  vtTarget.at(772).GalacticLat = 0.160689394015873;
  vtTarget.at(773).GalacticLat = -3.33451542532213;
  vtTarget.at(774).GalacticLat = 2.88245078980664;
  vtTarget.at(775).GalacticLat = -1.40175477099226;
  vtTarget.at(776).GalacticLat = 2.73491215909081;
  vtTarget.at(777).GalacticLat = 0.654721038205453;
  vtTarget.at(778).GalacticLat = 4.83599681160606;
  vtTarget.at(779).GalacticLat = 2.2416108504832;
  vtTarget.at(780).GalacticLat = -1.02471951290447;
  vtTarget.at(781).GalacticLat = -1.11370952241326;
  vtTarget.at(782).GalacticLat = 0.982714017762824;
  vtTarget.at(783).GalacticLat = 2.21113370592494;
  vtTarget.at(784).GalacticLat = 2.99000845222351;
  vtTarget.at(785).GalacticLat = -1.89441742959177;
  vtTarget.at(786).GalacticLat = -4.57851755366235;
  vtTarget.at(787).GalacticLat = 2.08996448278164;
  vtTarget.at(788).GalacticLat = -4.54186219761828;
  vtTarget.at(789).GalacticLat = 0.34904812028728;
  vtTarget.at(790).GalacticLat = -0.707001951886427;
  vtTarget.at(791).GalacticLat = -1.85785595746561;
  vtTarget.at(792).GalacticLat = -3.06741840286793;
  vtTarget.at(793).GalacticLat = -1.6090183157499;
  vtTarget.at(794).GalacticLat = 2.03417396597686;
  vtTarget.at(795).GalacticLat = 3.65071192995277;
  vtTarget.at(796).GalacticLat = -2.88166273350828;
  vtTarget.at(797).GalacticLat = 3.99615871299984;
  vtTarget.at(798).GalacticLat = -1.41523322792925;
  vtTarget.at(799).GalacticLat = -3.35600320710956;
  vtTarget.at(800).GalacticLat = -4.85696261155661;
  vtTarget.at(801).GalacticLat = -3.49541779872805;
  vtTarget.at(802).GalacticLat = -1.25639759333441;
  vtTarget.at(803).GalacticLat = -3.95158227098994;
  vtTarget.at(804).GalacticLat = 2.77643897122926;
  vtTarget.at(805).GalacticLat = -0.237090276031688;
  vtTarget.at(806).GalacticLat = 4.71071577294492;
  vtTarget.at(807).GalacticLat = 0.589200490612853;
  vtTarget.at(808).GalacticLat = -3.70115470294095;
  vtTarget.at(809).GalacticLat = 2.74455649481907;
  vtTarget.at(810).GalacticLat = -1.70962773592482;
  vtTarget.at(811).GalacticLat = 0.27849883116899;
  vtTarget.at(812).GalacticLat = -1.38302000325762;
  vtTarget.at(813).GalacticLat = -0.134937693450048;
  vtTarget.at(814).GalacticLat = -2.54647649825068;
  vtTarget.at(815).GalacticLat = -1.38189294340721;
  vtTarget.at(816).GalacticLat = -0.411344617463825;
  vtTarget.at(817).GalacticLat = 1.61696242562222;
  vtTarget.at(818).GalacticLat = -1.34801624999225;
  vtTarget.at(819).GalacticLat = -3.00457792515234;
  vtTarget.at(820).GalacticLat = 4.84558763604429;
  vtTarget.at(821).GalacticLat = 0.398354150552179;
  vtTarget.at(822).GalacticLat = 0.108587575916191;
  vtTarget.at(823).GalacticLat = -4.74521916786685;
  vtTarget.at(824).GalacticLat = 2.71350355902743;
  vtTarget.at(825).GalacticLat = 1.64472098262077;
  vtTarget.at(826).GalacticLat = -1.15152680901997;
  vtTarget.at(827).GalacticLat = 1.35810943920982;
  vtTarget.at(828).GalacticLat = 2.86165481190908;
  vtTarget.at(829).GalacticLat = -2.89989537126751;
  vtTarget.at(830).GalacticLat = 2.53933580193952;
  vtTarget.at(831).GalacticLat = 3.30669884525598;
  vtTarget.at(832).GalacticLat = -0.534639411992971;
  vtTarget.at(833).GalacticLat = -2.96028122268118;
  vtTarget.at(834).GalacticLat = 1.89628971606485;
  vtTarget.at(835).GalacticLat = -3.27601749523907;
  vtTarget.at(836).GalacticLat = 1.04706719592342;
  vtTarget.at(837).GalacticLat = -3.31224731542531;
  vtTarget.at(838).GalacticLat = 1.29060532725436;
  vtTarget.at(839).GalacticLat = -1.5250053943109;
  vtTarget.at(840).GalacticLat = 3.02441129347769;
  vtTarget.at(841).GalacticLat = -0.225631306109008;
  vtTarget.at(842).GalacticLat = 2.4790407844244;
  vtTarget.at(843).GalacticLat = -1.98364926550521;
  vtTarget.at(844).GalacticLat = 2.29258969621341;
  vtTarget.at(845).GalacticLat = 2.4433157022353;
  vtTarget.at(846).GalacticLat = -1.37314464132747;
  vtTarget.at(847).GalacticLat = 3.42125950073466;
  vtTarget.at(848).GalacticLat = -1.47209546456121;
  vtTarget.at(849).GalacticLat = -0.0534073307681737;
  vtTarget.at(850).GalacticLat = 4.3707169183708;
  vtTarget.at(851).GalacticLat = -2.62924926387602;
  vtTarget.at(852).GalacticLat = 3.30371893819824;
  vtTarget.at(853).GalacticLat = 2.89838229266758;
  vtTarget.at(854).GalacticLat = -2.07662089437395;
  vtTarget.at(855).GalacticLat = -3.72611231458552;
  vtTarget.at(856).GalacticLat = 0.869827365731609;
  vtTarget.at(857).GalacticLat = 3.05355809310398;
  vtTarget.at(858).GalacticLat = -0.666628698647096;
  vtTarget.at(859).GalacticLat = 4.70863634975807;
  vtTarget.at(860).GalacticLat = 0.0548113246394024;
  vtTarget.at(861).GalacticLat = -1.86316868327807;
  vtTarget.at(862).GalacticLat = -2.7237089199178;
  vtTarget.at(863).GalacticLat = -1.0019953305792;
  vtTarget.at(864).GalacticLat = 0.0125544774144029;
  vtTarget.at(865).GalacticLat = 3.09463219400438;
  vtTarget.at(866).GalacticLat = 0.0784766284756513;
  vtTarget.at(867).GalacticLat = -3.77149815544764;
  vtTarget.at(868).GalacticLat = -1.0947693097768;
  vtTarget.at(869).GalacticLat = -2.66421989548958;
  vtTarget.at(870).GalacticLat = -4.3810124607404;
  vtTarget.at(871).GalacticLat = -3.59856975207634;
  vtTarget.at(872).GalacticLat = 1.99056450258532;
  vtTarget.at(873).GalacticLat = 4.97671704542105;
  vtTarget.at(874).GalacticLat = 2.1359745468002;
  vtTarget.at(875).GalacticLat = 0.661269871257835;
  vtTarget.at(876).GalacticLat = 0.942862138807378;
  vtTarget.at(877).GalacticLat = 2.69991504286681;
  vtTarget.at(878).GalacticLat = -4.5764815499405;
  vtTarget.at(879).GalacticLat = -4.63868691000848;
  vtTarget.at(880).GalacticLat = 4.10926705308652;
  vtTarget.at(881).GalacticLat = -0.0896495414390328;
  vtTarget.at(882).GalacticLat = -4.7235141874338;
  vtTarget.at(883).GalacticLat = 0.477528482351989;
  vtTarget.at(884).GalacticLat = -0.757134738391176;
  vtTarget.at(885).GalacticLat = -2.80285140636484;
  vtTarget.at(886).GalacticLat = 0.239260520413513;
  vtTarget.at(887).GalacticLat = -3.58816380488016;
  vtTarget.at(888).GalacticLat = -1.6659394285556;
  vtTarget.at(889).GalacticLat = 2.59836614154225;
  vtTarget.at(890).GalacticLat = 3.74982985737081;
  vtTarget.at(891).GalacticLat = -4.18118409083733;
  vtTarget.at(892).GalacticLat = -3.62340463578662;
  vtTarget.at(893).GalacticLat = 3.09091314692303;
  vtTarget.at(894).GalacticLat = 0.937441560277722;
  vtTarget.at(895).GalacticLat = 4.73340527989594;
  vtTarget.at(896).GalacticLat = -2.40832451890167;
  vtTarget.at(897).GalacticLat = -3.58113131307672;
  vtTarget.at(898).GalacticLat = -4.85353037180721;
  vtTarget.at(899).GalacticLat = -0.808557078302673;
  vtTarget.at(900).GalacticLat = -4.39455371886985;
  vtTarget.at(901).GalacticLat = -1.87914125909114;
  vtTarget.at(902).GalacticLat = -2.03884339692896;
  vtTarget.at(903).GalacticLat = 4.98606635454499;
  vtTarget.at(904).GalacticLat = -1.38757635151006;
  vtTarget.at(905).GalacticLat = -4.68316521978571;
  vtTarget.at(906).GalacticLat = 0.171930420537912;
  vtTarget.at(907).GalacticLat = 3.18042796451645;
  vtTarget.at(908).GalacticLat = -0.522949706260287;
  vtTarget.at(909).GalacticLat = 3.82424642026295;
  vtTarget.at(910).GalacticLat = -2.14724014825678;
  vtTarget.at(911).GalacticLat = 0.888921981988631;
  vtTarget.at(912).GalacticLat = -1.48219305308569;
  vtTarget.at(913).GalacticLat = 4.46242320117486;
  vtTarget.at(914).GalacticLat = 1.28630399851544;
  vtTarget.at(915).GalacticLat = 0.0743908974772332;
  vtTarget.at(916).GalacticLat = -2.0949352808638;
  vtTarget.at(917).GalacticLat = 4.22109043559383;
  vtTarget.at(918).GalacticLat = -3.5386181047872;
  vtTarget.at(919).GalacticLat = 3.03214204583376;
  vtTarget.at(920).GalacticLat = -2.20049637350416;
  vtTarget.at(921).GalacticLat = -2.03486576077072;
  vtTarget.at(922).GalacticLat = 3.77280946376902;
  vtTarget.at(923).GalacticLat = -1.69385926763292;
  vtTarget.at(924).GalacticLat = 2.46777784042451;
  vtTarget.at(925).GalacticLat = -2.43964912438696;
  vtTarget.at(926).GalacticLat = -2.67438078782555;
  vtTarget.at(927).GalacticLat = -3.28475760976371;
  vtTarget.at(928).GalacticLat = -4.96777403004704;
  vtTarget.at(929).GalacticLat = -4.81439704895471;
  vtTarget.at(930).GalacticLat = -0.521264252784446;
  vtTarget.at(931).GalacticLat = 0.0290511873260062;
  vtTarget.at(932).GalacticLat = -1.18901141447104;
  vtTarget.at(933).GalacticLat = 3.71743917454181;
  vtTarget.at(934).GalacticLat = 2.10760518486132;
  vtTarget.at(935).GalacticLat = -4.87362615235055;
  vtTarget.at(936).GalacticLat = 4.11933962761915;
  vtTarget.at(937).GalacticLat = 3.31000938626243;
  vtTarget.at(938).GalacticLat = 3.80962568005141;
  vtTarget.at(939).GalacticLat = -3.17088225679021;
  vtTarget.at(940).GalacticLat = 3.7893306333709;
  vtTarget.at(941).GalacticLat = -0.561240586765806;
  vtTarget.at(942).GalacticLat = -3.6968394459098;
  vtTarget.at(943).GalacticLat = -2.70154797920243;
  vtTarget.at(944).GalacticLat = 1.68137280814854;
  vtTarget.at(945).GalacticLat = -0.506702009587737;
  vtTarget.at(946).GalacticLat = -4.5449140648098;
  vtTarget.at(947).GalacticLat = 1.67765695072147;
  vtTarget.at(948).GalacticLat = 2.99271674932524;
  vtTarget.at(949).GalacticLat = 3.30322605385514;
  vtTarget.at(950).GalacticLat = 3.8078511536762;
  vtTarget.at(951).GalacticLat = -1.73984506488872;
  vtTarget.at(952).GalacticLat = 3.37692053915579;
  vtTarget.at(953).GalacticLat = -1.11162287391401;
  vtTarget.at(954).GalacticLat = -3.13800306531155;
  vtTarget.at(955).GalacticLat = 0.380307577077045;
  vtTarget.at(956).GalacticLat = 2.52964861892595;
  vtTarget.at(957).GalacticLat = 4.40727407726562;
  vtTarget.at(958).GalacticLat = 3.26387762002186;
  vtTarget.at(959).GalacticLat = 4.29244752054712;
  vtTarget.at(960).GalacticLat = -4.75144001257462;
  vtTarget.at(961).GalacticLat = 3.13197842194512;
  vtTarget.at(962).GalacticLat = -4.6451134175622;
  vtTarget.at(963).GalacticLat = -2.43730479327713;
  vtTarget.at(964).GalacticLat = -1.7676420068279;
  vtTarget.at(965).GalacticLat = 0.282537582939334;
  vtTarget.at(966).GalacticLat = 2.25124859340088;
  vtTarget.at(967).GalacticLat = 3.12883268488705;
  vtTarget.at(968).GalacticLat = -4.10302868684689;
  vtTarget.at(969).GalacticLat = -1.82689063098862;
  vtTarget.at(970).GalacticLat = 3.57978206350211;
  vtTarget.at(971).GalacticLat = 1.81673826676646;
  vtTarget.at(972).GalacticLat = -1.10637469448301;
  vtTarget.at(973).GalacticLat = 1.41039297050518;
  vtTarget.at(974).GalacticLat = 2.78471921225995;
  vtTarget.at(975).GalacticLat = -2.23868710226379;
  vtTarget.at(976).GalacticLat = -2.07245490372278;
  vtTarget.at(977).GalacticLat = 2.43133965670218;
  vtTarget.at(978).GalacticLat = 2.49403997520638;
  vtTarget.at(979).GalacticLat = 2.19678800208721;
  vtTarget.at(980).GalacticLat = 3.107712895825;
  vtTarget.at(981).GalacticLat = 1.27022542752512;
  vtTarget.at(982).GalacticLat = 2.73035372955047;
  vtTarget.at(983).GalacticLat = -1.07664309750931;
  vtTarget.at(984).GalacticLat = 2.02308943371608;
  vtTarget.at(985).GalacticLat = -1.87078981246401;
  vtTarget.at(986).GalacticLat = -1.09762429769335;
  vtTarget.at(987).GalacticLat = -4.91851913946133;
  vtTarget.at(988).GalacticLat = -0.367194285541506;
  vtTarget.at(989).GalacticLat = 1.11866896780575;
  vtTarget.at(990).GalacticLat = 4.08706573344312;
  vtTarget.at(991).GalacticLat = -0.357098967773751;
  vtTarget.at(992).GalacticLat = -1.9814774749148;
  vtTarget.at(993).GalacticLat = 2.8986847324737;
  vtTarget.at(994).GalacticLat = 0.51162048609235;
  vtTarget.at(995).GalacticLat = 3.64694538128943;
  vtTarget.at(996).GalacticLat = -0.0140025406881339;
  vtTarget.at(997).GalacticLat = -3.56765224757409;
  vtTarget.at(998).GalacticLat = 3.26082091921149;
  vtTarget.at(999).GalacticLat = -3.04358727872877;
  vtTarget.at(1000).GalacticLat = 0.551332891056253;
  vtTarget.at(1001).GalacticLat = -4.87504932392918;
  vtTarget.at(1002).GalacticLat = 0.216378831724589;
  vtTarget.at(1003).GalacticLat = 4.86747062054854;
  vtTarget.at(1004).GalacticLat = 4.2399591683611;
  vtTarget.at(1005).GalacticLat = -4.94070729377277;
  vtTarget.at(1006).GalacticLat = 1.1831280202489;
  vtTarget.at(1007).GalacticLat = 4.93209547171226;
  vtTarget.at(1008).GalacticLat = -2.0558532045382;
  vtTarget.at(1009).GalacticLat = 4.24290800718273;
  vtTarget.at(1010).GalacticLat = -4.09988082818232;
  vtTarget.at(1011).GalacticLat = 1.35659099194529;
  vtTarget.at(1012).GalacticLat = -3.50994242185727;
  vtTarget.at(1013).GalacticLat = 3.86015775463882;
  vtTarget.at(1014).GalacticLat = -3.74309480789118;
  vtTarget.at(1015).GalacticLat = -2.29423006078898;
  vtTarget.at(1016).GalacticLat = -4.09260507439059;
  vtTarget.at(1017).GalacticLat = 2.93781281384845;
  vtTarget.at(1018).GalacticLat = -0.584686693982919;
  vtTarget.at(1019).GalacticLat = 0.937782133962768;
  vtTarget.at(1020).GalacticLat = -0.279856132690075;
  vtTarget.at(1021).GalacticLat = -1.59890118644358;
  vtTarget.at(1022).GalacticLat = -1.12055270448456;
  vtTarget.at(1023).GalacticLat = -2.78801631666033;
  vtTarget.at(1024).GalacticLat = -3.18014824515097;
  vtTarget.at(1025).GalacticLat = 2.46391789750996;
  vtTarget.at(1026).GalacticLat = 1.71311980549782;
  vtTarget.at(1027).GalacticLat = 3.74648233353032;
  vtTarget.at(1028).GalacticLat = 2.08927937290351;
  vtTarget.at(1029).GalacticLat = -2.27970582332532;
  vtTarget.at(1030).GalacticLat = 2.92805632552746;
  vtTarget.at(1031).GalacticLat = -4.39791313561392;
  vtTarget.at(1032).GalacticLat = 4.94587296210509;
  vtTarget.at(1033).GalacticLat = 0.679978550401081;
  vtTarget.at(1034).GalacticLat = -2.48501015914394;
  vtTarget.at(1035).GalacticLat = 0.0225210125212785;
  vtTarget.at(1036).GalacticLat = 0.542993997272596;
  vtTarget.at(1037).GalacticLat = -0.430204991045335;
  vtTarget.at(1038).GalacticLat = -2.61948570931159;
  vtTarget.at(1039).GalacticLat = -2.8013692701353;
  vtTarget.at(1040).GalacticLat = -0.48501602373308;
  vtTarget.at(1041).GalacticLat = -0.762429117861779;
  vtTarget.at(1042).GalacticLat = -4.37894553704797;
  vtTarget.at(1043).GalacticLat = 4.2642976168921;
  vtTarget.at(1044).GalacticLat = 1.85435976017722;
  vtTarget.at(1045).GalacticLat = -1.49452010101944;
  vtTarget.at(1046).GalacticLat = 0.580820748345047;
  vtTarget.at(1047).GalacticLat = 2.74694788879915;
  vtTarget.at(1048).GalacticLat = -0.855920959906338;
  vtTarget.at(1049).GalacticLat = -1.02528562219566;
  vtTarget.at(1050).GalacticLat = -1.70630215154866;
  vtTarget.at(1051).GalacticLat = -3.95079629817555;
  vtTarget.at(1052).GalacticLat = 1.25731039045941;
  vtTarget.at(1053).GalacticLat = 1.84972295713706;
  vtTarget.at(1054).GalacticLat = 2.54847123650532;
  vtTarget.at(1055).GalacticLat = 0.40261975607662;
  vtTarget.at(1056).GalacticLat = -2.70693835391513;
  vtTarget.at(1057).GalacticLat = -2.94578725812502;
  vtTarget.at(1058).GalacticLat = 0.2286167244755;
  vtTarget.at(1059).GalacticLat = -1.34027585967168;
  vtTarget.at(1060).GalacticLat = 1.60993767318648;
  vtTarget.at(1061).GalacticLat = -4.39543496365338;
  vtTarget.at(1062).GalacticLat = -3.81970714102049;
  vtTarget.at(1063).GalacticLat = 1.23239898322596;
  vtTarget.at(1064).GalacticLat = -0.591515343271026;
  vtTarget.at(1065).GalacticLat = 3.50036826850956;
  vtTarget.at(1066).GalacticLat = 1.048099725173;
  vtTarget.at(1067).GalacticLat = 3.31069552260136;
  vtTarget.at(1068).GalacticLat = -3.08479778238602;
  vtTarget.at(1069).GalacticLat = 1.71087958010294;
  vtTarget.at(1070).GalacticLat = -1.65372196704207;
  vtTarget.at(1071).GalacticLat = -0.84245106954728;
  vtTarget.at(1072).GalacticLat = -4.29471549088293;
  vtTarget.at(1073).GalacticLat = 0.739395812349084;
  vtTarget.at(1074).GalacticLat = 2.22969601540033;
  vtTarget.at(1075).GalacticLat = 1.7672284556264;
  vtTarget.at(1076).GalacticLat = 2.86284167937682;
  vtTarget.at(1077).GalacticLat = -4.69519971560762;
  vtTarget.at(1078).GalacticLat = -3.78676121729523;
  vtTarget.at(1079).GalacticLat = -3.11454241391354;
  vtTarget.at(1080).GalacticLat = 2.54381175027653;
  vtTarget.at(1081).GalacticLat = -2.85142860892122;
  vtTarget.at(1082).GalacticLat = 1.25477427700582;
  vtTarget.at(1083).GalacticLat = 0.626194461932465;
  vtTarget.at(1084).GalacticLat = -0.993725336487739;
  vtTarget.at(1085).GalacticLat = 1.29037503526203;
  vtTarget.at(1086).GalacticLat = -4.4420460041451;
  vtTarget.at(1087).GalacticLat = -1.42202066992838;
  vtTarget.at(1088).GalacticLat = 4.82015103056132;
  vtTarget.at(1089).GalacticLat = -2.47139974945459;
  vtTarget.at(1090).GalacticLat = -2.53342712758706;
  vtTarget.at(1091).GalacticLat = 3.96691710729028;
  vtTarget.at(1092).GalacticLat = -3.89222478463451;
  vtTarget.at(1093).GalacticLat = -0.495739103866711;
  vtTarget.at(1094).GalacticLat = -1.83122840364841;
  vtTarget.at(1095).GalacticLat = 3.10861254304131;
  vtTarget.at(1096).GalacticLat = -2.62045834108515;
  vtTarget.at(1097).GalacticLat = -3.02688549971056;
  vtTarget.at(1098).GalacticLat = -3.51324977384928;
  vtTarget.at(1099).GalacticLat = -1.27459582896514;

  vtTarget.at(0).GalacticLong = -0.571221976398421;
  vtTarget.at(1).GalacticLong = -0.171056020130747;
  vtTarget.at(2).GalacticLong = -0.552686695343407;
  vtTarget.at(3).GalacticLong = 0.584813710871605;
  vtTarget.at(4).GalacticLong = 0.963373062814946;
  vtTarget.at(5).GalacticLong = -0.340708562551168;
  vtTarget.at(6).GalacticLong = -0.543118120378066;
  vtTarget.at(7).GalacticLong = -0.700864772758388;
  vtTarget.at(8).GalacticLong = 0.889843768545944;
  vtTarget.at(9).GalacticLong = 0.00869345457624546;
  vtTarget.at(10).GalacticLong = -0.0860552720687565;
  vtTarget.at(11).GalacticLong = 0.0516106168165547;
  vtTarget.at(12).GalacticLong = 0.0508229573356329;
  vtTarget.at(13).GalacticLong = -0.299176252068633;
  vtTarget.at(14).GalacticLong = 0.292185814324491;
  vtTarget.at(15).GalacticLong = 1.02246045029859;
  vtTarget.at(16).GalacticLong = 0.586409853548903;
  vtTarget.at(17).GalacticLong = -0.0106618225528939;
  vtTarget.at(18).GalacticLong = -0.0348674821528791;
  vtTarget.at(19).GalacticLong = -0.544672980609567;
  vtTarget.at(20).GalacticLong = -1.31835811802153;
  vtTarget.at(21).GalacticLong = -0.812893645458297;
  vtTarget.at(22).GalacticLong = -0.836429043036418;
  vtTarget.at(23).GalacticLong = -1.44847593950578;
  vtTarget.at(24).GalacticLong = -0.945295531267032;
  vtTarget.at(25).GalacticLong = -2.04656701030529;
  vtTarget.at(26).GalacticLong = -1.19191606496102;
  vtTarget.at(27).GalacticLong = -0.730367146451408;
  vtTarget.at(28).GalacticLong = -1.55355392805346;
  vtTarget.at(29).GalacticLong = -1.06899522112577;
  vtTarget.at(30).GalacticLong = -0.355212856125503;
  vtTarget.at(31).GalacticLong = -0.480244053551432;
  vtTarget.at(32).GalacticLong = -1.06608409101584;
  vtTarget.at(33).GalacticLong = -0.171333758437633;
  vtTarget.at(34).GalacticLong = -0.290495854088603;
  vtTarget.at(35).GalacticLong = -0.512020541948482;
  vtTarget.at(36).GalacticLong = -1.01149709105066;
  vtTarget.at(37).GalacticLong = 3.30090183790877;
  vtTarget.at(38).GalacticLong = 3.52863771982634;
  vtTarget.at(39).GalacticLong = 2.55775129871701;
  vtTarget.at(40).GalacticLong = 2.55954593506835;
  vtTarget.at(41).GalacticLong = 4.02497475202947;
  vtTarget.at(42).GalacticLong = 2.76694581655972;
  vtTarget.at(43).GalacticLong = 2.67074069390283;
  vtTarget.at(44).GalacticLong = 2.90394611180939;
  vtTarget.at(45).GalacticLong = 4.12449582770493;
  vtTarget.at(46).GalacticLong = 4.1851268985903;
  vtTarget.at(47).GalacticLong = 3.33278797235773;
  vtTarget.at(48).GalacticLong = 2.7554560711728;
  vtTarget.at(49).GalacticLong = -0.353416783840857;
  vtTarget.at(50).GalacticLong = -0.284735660041338;
  vtTarget.at(51).GalacticLong = -0.265830490803246;
  vtTarget.at(52).GalacticLong = -0.325048094019605;
  vtTarget.at(53).GalacticLong = -0.231675969033742;
  vtTarget.at(54).GalacticLong = -0.241501558870553;
  vtTarget.at(55).GalacticLong = -0.382473162777631;
  vtTarget.at(56).GalacticLong = 0.00217660801164901;
  vtTarget.at(57).GalacticLong = -0.489563089929781;
  vtTarget.at(58).GalacticLong = 0.0220469825851477;
  vtTarget.at(59).GalacticLong = -0.857408051736337;
  vtTarget.at(60).GalacticLong = -0.631795316876644;
  vtTarget.at(61).GalacticLong = -4.48053626120991;
  vtTarget.at(62).GalacticLong = -4.97758824571904;
  vtTarget.at(63).GalacticLong = -4.10527835575641;
  vtTarget.at(64).GalacticLong = -4.68377181597608;
  vtTarget.at(65).GalacticLong = -4.79579315992014;
  vtTarget.at(66).GalacticLong = -4.7816868404566;
  vtTarget.at(67).GalacticLong = -4.4436085020193;
  vtTarget.at(68).GalacticLong = -3.8264480040255;
  vtTarget.at(69).GalacticLong = -4.03766547661004;
  vtTarget.at(70).GalacticLong = -4.16451200017516;
  vtTarget.at(71).GalacticLong = -3.76519724645186;
  vtTarget.at(72).GalacticLong = -4.00335328817032;
  vtTarget.at(73).GalacticLong = -2.67201205747468;
  vtTarget.at(74).GalacticLong = -2.27720500837858;
  vtTarget.at(75).GalacticLong = -2.31444243343463;
  vtTarget.at(76).GalacticLong = -2.30821603331922;
  vtTarget.at(77).GalacticLong = -2.23139626059439;
  vtTarget.at(78).GalacticLong = -3.86294579593899;
  vtTarget.at(79).GalacticLong = -3.28539738671735;
  vtTarget.at(80).GalacticLong = -3.59645861412712;
  vtTarget.at(81).GalacticLong = -4.08265725203623;
  vtTarget.at(82).GalacticLong = -4.198954588176;
  vtTarget.at(83).GalacticLong = -3.83430167488198;
  vtTarget.at(84).GalacticLong = -4.01230543608727;
  vtTarget.at(85).GalacticLong = -4.06825323737586;
  vtTarget.at(86).GalacticLong = -3.77978199796322;
  vtTarget.at(87).GalacticLong = -4.78930139397619;
  vtTarget.at(88).GalacticLong = -4.52969314634881;
  vtTarget.at(89).GalacticLong = -4.69036777064734;
  vtTarget.at(90).GalacticLong = -4.73797280324936;
  vtTarget.at(91).GalacticLong = -4.86638787434943;
  vtTarget.at(92).GalacticLong = -3.37574399595229;
  vtTarget.at(93).GalacticLong = -3.52591952972716;
  vtTarget.at(94).GalacticLong = -4.74888765895234;
  vtTarget.at(95).GalacticLong = -3.6883054134535;
  vtTarget.at(96).GalacticLong = -5.01219219129742;
  vtTarget.at(97).GalacticLong = -4.02055665491743;
  vtTarget.at(98).GalacticLong = -4.26574336556421;
  vtTarget.at(99).GalacticLong = -4.62038316251741;
  vtTarget.at(100).GalacticLong = -4.39396749232844;
  vtTarget.at(101).GalacticLong = -4.67670225826798;
  vtTarget.at(102).GalacticLong = -4.71292835160637;
  vtTarget.at(103).GalacticLong = -4.740711086213;
  vtTarget.at(104).GalacticLong = -3.42837807856279;
  vtTarget.at(105).GalacticLong = -4.39900667002803;
  vtTarget.at(106).GalacticLong = -3.61261780447219;
  vtTarget.at(107).GalacticLong = -4.05530329873478;
  vtTarget.at(108).GalacticLong = -3.5211884067224;
  vtTarget.at(109).GalacticLong = -3.5167696831986;
  vtTarget.at(110).GalacticLong = -3.88780526171506;
  vtTarget.at(111).GalacticLong = 2.4374291036348;
  vtTarget.at(112).GalacticLong = 1.8409622777816;
  vtTarget.at(113).GalacticLong = 2.58196185947658;
  vtTarget.at(114).GalacticLong = 1.92666749407604;
  vtTarget.at(115).GalacticLong = 2.74852497986536;
  vtTarget.at(116).GalacticLong = 2.54829570452276;
  vtTarget.at(117).GalacticLong = 3.24054698315542;
  vtTarget.at(118).GalacticLong = 2.21261318462785;
  vtTarget.at(119).GalacticLong = 2.62753976887988;
  vtTarget.at(120).GalacticLong = 2.58962234774736;
  vtTarget.at(121).GalacticLong = 2.50171241179295;
  vtTarget.at(122).GalacticLong = 2.43342803777289;
  vtTarget.at(123).GalacticLong = 1.36881669719178;
  vtTarget.at(124).GalacticLong = 1.36985455275222;
  vtTarget.at(125).GalacticLong = 2.85044342443211;
  vtTarget.at(126).GalacticLong = 2.71299242764163;
  vtTarget.at(127).GalacticLong = 2.4719834721672;
  vtTarget.at(128).GalacticLong = 2.48338229835854;
  vtTarget.at(129).GalacticLong = 2.29612619797518;
  vtTarget.at(130).GalacticLong = 1.82902319311853;
  vtTarget.at(131).GalacticLong = 1.71549835452706;
  vtTarget.at(132).GalacticLong = 1.63935528202787;
  vtTarget.at(133).GalacticLong = 1.58894344564995;
  vtTarget.at(134).GalacticLong = 1.41383499127863;
  vtTarget.at(135).GalacticLong = 1.46077285314476;
  vtTarget.at(136).GalacticLong = 3.03670474895441;
  vtTarget.at(137).GalacticLong = 2.99355085050277;
  vtTarget.at(138).GalacticLong = 2.40796196321659;
  vtTarget.at(139).GalacticLong = 2.20138902818764;
  vtTarget.at(140).GalacticLong = 2.65497684462471;
  vtTarget.at(141).GalacticLong = 2.45830262497527;
  vtTarget.at(142).GalacticLong = 2.11228770706201;
  vtTarget.at(143).GalacticLong = 1.77476348505983;
  vtTarget.at(144).GalacticLong = 1.76290918454467;
  vtTarget.at(145).GalacticLong = 2.68049019386842;
  vtTarget.at(146).GalacticLong = 1.90980214401138;
  vtTarget.at(147).GalacticLong = 1.92267743145946;
  vtTarget.at(148).GalacticLong = 1.82649475442465;
  vtTarget.at(149).GalacticLong = 2.7306408842349;
  vtTarget.at(150).GalacticLong = 2.59119178301227;
  vtTarget.at(151).GalacticLong = 2.60099888058925;
  vtTarget.at(152).GalacticLong = 2.60914408733404;
  vtTarget.at(153).GalacticLong = 2.35707021015541;
  vtTarget.at(154).GalacticLong = 2.27734604731215;
  vtTarget.at(155).GalacticLong = 2.23412473950487;
  vtTarget.at(156).GalacticLong = 1.81928311307979;
  vtTarget.at(157).GalacticLong = 1.7478597580552;
  vtTarget.at(158).GalacticLong = 2.48140122398399;
  vtTarget.at(159).GalacticLong = 2.60592734026283;
  vtTarget.at(160).GalacticLong = 3.02403910983131;
  vtTarget.at(161).GalacticLong = 1.65982155917886;
  vtTarget.at(162).GalacticLong = 2.12412665948226;
  vtTarget.at(163).GalacticLong = 1.93447490594772;
  vtTarget.at(164).GalacticLong = 1.9530325323211;
  vtTarget.at(165).GalacticLong = 1.58221505595911;
  vtTarget.at(166).GalacticLong = 2.38595359956626;
  vtTarget.at(167).GalacticLong = 2.28394173352815;
  vtTarget.at(168).GalacticLong = 2.34610587175759;
  vtTarget.at(169).GalacticLong = 2.28105990828908;
  vtTarget.at(170).GalacticLong = 2.05705387039846;
  vtTarget.at(171).GalacticLong = 1.60522052131146;
  vtTarget.at(172).GalacticLong = 2.67339907532209;
  vtTarget.at(173).GalacticLong = 2.61048929855139;
  vtTarget.at(174).GalacticLong = 2.43809582912083;
  vtTarget.at(175).GalacticLong = 2.24558382541602;
  vtTarget.at(176).GalacticLong = 2.13810435054878;
  vtTarget.at(177).GalacticLong = 2.8170554766582;
  vtTarget.at(178).GalacticLong = 2.28423017792823;
  vtTarget.at(179).GalacticLong = 1.83919676736202;
  vtTarget.at(180).GalacticLong = 2.94584667975335;
  vtTarget.at(181).GalacticLong = 3.93008140290289;
  vtTarget.at(182).GalacticLong = 4.54582250126185;
  vtTarget.at(183).GalacticLong = 4.13765888044207;
  vtTarget.at(184).GalacticLong = 4.20404485481185;
  vtTarget.at(185).GalacticLong = 4.06088059227523;
  vtTarget.at(186).GalacticLong = 4.55150943189379;
  vtTarget.at(187).GalacticLong = 4.42210288448753;
  vtTarget.at(188).GalacticLong = 4.25089458590862;
  vtTarget.at(189).GalacticLong = 4.15352881774284;
  vtTarget.at(190).GalacticLong = -0.9543309034414;
  vtTarget.at(191).GalacticLong = -0.151082312682952;
  vtTarget.at(192).GalacticLong = 0.559954392144918;
  vtTarget.at(193).GalacticLong = 0.220711131153594;
  vtTarget.at(194).GalacticLong = 0.981406207423325;
  vtTarget.at(195).GalacticLong = 0.89475025360862;
  vtTarget.at(196).GalacticLong = -2.52090605532933;
  vtTarget.at(197).GalacticLong = -2.36034056809797;
  vtTarget.at(198).GalacticLong = -2.0273358118864;
  vtTarget.at(199).GalacticLong = -2.01140230318011;
  vtTarget.at(200).GalacticLong = -2.59671796129189;
  vtTarget.at(201).GalacticLong = -2.69156264971134;
  vtTarget.at(202).GalacticLong = -2.81830202756116;
  vtTarget.at(203).GalacticLong = -4.88242753672233;
  vtTarget.at(204).GalacticLong = -4.34245831073652;
  vtTarget.at(205).GalacticLong = -4.2633065542896;
  vtTarget.at(206).GalacticLong = 4.88695235030289;
  vtTarget.at(207).GalacticLong = 4.07766984305543;
  vtTarget.at(208).GalacticLong = 2.0941623840985;
  vtTarget.at(209).GalacticLong = 2.39370406432256;
  vtTarget.at(210).GalacticLong = 0.605144058655807;
  vtTarget.at(211).GalacticLong = 0.445271934449683;
  vtTarget.at(212).GalacticLong = 0.805941486063945;
  vtTarget.at(213).GalacticLong = 0.341970602421905;
  vtTarget.at(214).GalacticLong = 2.17833473789806;
  vtTarget.at(215).GalacticLong = -4.46628097939342;
  vtTarget.at(216).GalacticLong = -2.28985739539738;
  vtTarget.at(217).GalacticLong = -4.59172298450752;
  vtTarget.at(218).GalacticLong = 0.583309360686938;
  vtTarget.at(219).GalacticLong = 0.541330049152442;
  vtTarget.at(220).GalacticLong = 0.409797723314372;
  vtTarget.at(221).GalacticLong = 0.129547209919364;
  vtTarget.at(222).GalacticLong = -3.21692892235877;
  vtTarget.at(223).GalacticLong = -3.42149193613681;
  vtTarget.at(224).GalacticLong = -3.81228571076501;
  vtTarget.at(225).GalacticLong = 1.29869941729723;
  vtTarget.at(226).GalacticLong = 0.595470689800515;
  vtTarget.at(227).GalacticLong = -3.15225136753758;
  vtTarget.at(228).GalacticLong = -4.15137809916473;
  vtTarget.at(229).GalacticLong = -4.22176027507017;
  vtTarget.at(230).GalacticLong = 1.43763836129132;
  vtTarget.at(231).GalacticLong = 1.4258495148651;
  vtTarget.at(232).GalacticLong = 0.829106662501559;
  vtTarget.at(233).GalacticLong = 0.867824129975193;
  vtTarget.at(234).GalacticLong = 0.63093768227031;
  vtTarget.at(235).GalacticLong = 0.346941271694604;
  vtTarget.at(236).GalacticLong = -2.13532543175109;
  vtTarget.at(237).GalacticLong = 1.38939460837827;
  vtTarget.at(238).GalacticLong = 0.956202418877378;
  vtTarget.at(239).GalacticLong = 0.856667748861638;
  vtTarget.at(240).GalacticLong = 0.89239650916279;
  vtTarget.at(241).GalacticLong = 0.598552215586833;
  vtTarget.at(242).GalacticLong = -0.0216455059095955;
  vtTarget.at(243).GalacticLong = 1.37983285227545;
  vtTarget.at(244).GalacticLong = 0.738560389055428;
  vtTarget.at(245).GalacticLong = 1.15713636781977;
  vtTarget.at(246).GalacticLong = 1.01138722581052;
  vtTarget.at(247).GalacticLong = 1.31997213858924;
  vtTarget.at(248).GalacticLong = 1.2047393613269;
  vtTarget.at(249).GalacticLong = 0.708075899787784;
  vtTarget.at(250).GalacticLong = -0.55601746497572;
  vtTarget.at(251).GalacticLong = -0.57397216846817;
  vtTarget.at(252).GalacticLong = 1.23748359722623;
  vtTarget.at(253).GalacticLong = 1.18716543026798;
  vtTarget.at(254).GalacticLong = 1.07949724747705;
  vtTarget.at(255).GalacticLong = 0.989730241065897;
  vtTarget.at(256).GalacticLong = 1.47234960379022;
  vtTarget.at(257).GalacticLong = 1.25688307717995;
  vtTarget.at(258).GalacticLong = 1.17839153160888;
  vtTarget.at(259).GalacticLong = 1.59917471468481;
  vtTarget.at(260).GalacticLong = 1.54392738534772;
  vtTarget.at(261).GalacticLong = -1.17980113621982;
  vtTarget.at(262).GalacticLong = -1.17722596213394;
  vtTarget.at(263).GalacticLong = 1.91261736173068;
  vtTarget.at(264).GalacticLong = 1.88599566348811;
  vtTarget.at(265).GalacticLong = 1.88456609059748;
  vtTarget.at(266).GalacticLong = 3.55797381745236;
  vtTarget.at(267).GalacticLong = 3.55407041716134;
  vtTarget.at(268).GalacticLong = 2.19689577130936;
  vtTarget.at(269).GalacticLong = 2.13157619665938;
  vtTarget.at(270).GalacticLong = 2.57580079697804;
  vtTarget.at(271).GalacticLong = 2.50020792959624;
  vtTarget.at(272).GalacticLong = 4.00876829467294;
  vtTarget.at(273).GalacticLong = 3.16053534741641;
  vtTarget.at(274).GalacticLong = 3.03753455958397;
  vtTarget.at(275).GalacticLong = 3.05496854436014;
  vtTarget.at(276).GalacticLong = 3.6251655430263;
  vtTarget.at(277).GalacticLong = -0.571418729344014;
  vtTarget.at(278).GalacticLong = -1.10523210539767;
  vtTarget.at(279).GalacticLong = -2.42572783796597;
  vtTarget.at(280).GalacticLong = 0.527421180134272;
  vtTarget.at(281).GalacticLong = 0.720649759411186;
  vtTarget.at(282).GalacticLong = -0.248105935028483;
  vtTarget.at(283).GalacticLong = -3.99313086056283;
  vtTarget.at(284).GalacticLong = 0.511632643536916;
  vtTarget.at(285).GalacticLong = -4.7804497114889;
  vtTarget.at(286).GalacticLong = -2.54703913247528;
  vtTarget.at(287).GalacticLong = -1.22035449601156;
  vtTarget.at(288).GalacticLong = 2.7053747463968;
  vtTarget.at(289).GalacticLong = -2.71360509803532;
  vtTarget.at(290).GalacticLong = -1.2169790744275;
  vtTarget.at(291).GalacticLong = 3.12703259025102;
  vtTarget.at(292).GalacticLong = 3.60891124029475;
  vtTarget.at(293).GalacticLong = 0.883475269170821;
  vtTarget.at(294).GalacticLong = 0.0578585991925697;
  vtTarget.at(295).GalacticLong = -3.51993944123976;
  vtTarget.at(296).GalacticLong = -0.667378524632819;
  vtTarget.at(297).GalacticLong = -1.18125214565619;
  vtTarget.at(298).GalacticLong = -4.13208426531788;
  vtTarget.at(299).GalacticLong = 4.47000479276638;
  vtTarget.at(300).GalacticLong = 0.54761980262186;
  vtTarget.at(301).GalacticLong = -3.5965644393252;
  vtTarget.at(302).GalacticLong = 2.19106658881847;
  vtTarget.at(303).GalacticLong = -0.786610138879369;
  vtTarget.at(304).GalacticLong = -4.75595341472029;
  vtTarget.at(305).GalacticLong = 4.39156072797074;
  vtTarget.at(306).GalacticLong = -4.89176257249132;
  vtTarget.at(307).GalacticLong = -4.91076910137818;
  vtTarget.at(308).GalacticLong = -2.4312907107254;
  vtTarget.at(309).GalacticLong = 2.35676454211655;
  vtTarget.at(310).GalacticLong = -3.10619575744341;
  vtTarget.at(311).GalacticLong = 1.29725066320515;
  vtTarget.at(312).GalacticLong = -1.05026004207313;
  vtTarget.at(313).GalacticLong = 1.32121416547402;
  vtTarget.at(314).GalacticLong = -3.68643504623033;
  vtTarget.at(315).GalacticLong = 5.44065413112804E-05;
  vtTarget.at(316).GalacticLong = -2.65625789672429;
  vtTarget.at(317).GalacticLong = -4.80235203539746;
  vtTarget.at(318).GalacticLong = 1.36927844969626;
  vtTarget.at(319).GalacticLong = -3.05788110029467;
  vtTarget.at(320).GalacticLong = 4.8635703749982;
  vtTarget.at(321).GalacticLong = -2.35627733473898;
  vtTarget.at(322).GalacticLong = 1.36699150865568;
  vtTarget.at(323).GalacticLong = 4.56482190312673;
  vtTarget.at(324).GalacticLong = 0.590202336164082;
  vtTarget.at(325).GalacticLong = 2.82759608603886;
  vtTarget.at(326).GalacticLong = -0.183359479694877;
  vtTarget.at(327).GalacticLong = -3.54989600893066;
  vtTarget.at(328).GalacticLong = 1.08455938247158;
  vtTarget.at(329).GalacticLong = -2.68323955485389;
  vtTarget.at(330).GalacticLong = 4.16972111738318;
  vtTarget.at(331).GalacticLong = 0.619680766973817;
  vtTarget.at(332).GalacticLong = 4.6033610871961;
  vtTarget.at(333).GalacticLong = 1.16586601896645;
  vtTarget.at(334).GalacticLong = -2.10098010311726;
  vtTarget.at(335).GalacticLong = 2.47712178838168;
  vtTarget.at(336).GalacticLong = -4.9795290399367;
  vtTarget.at(337).GalacticLong = -3.15879452185361;
  vtTarget.at(338).GalacticLong = 4.11258070019275;
  vtTarget.at(339).GalacticLong = -1.61104899291426;
  vtTarget.at(340).GalacticLong = -3.78102569334462;
  vtTarget.at(341).GalacticLong = 1.02943191878362;
  vtTarget.at(342).GalacticLong = 1.61484993202823;
  vtTarget.at(343).GalacticLong = 3.29611596453826;
  vtTarget.at(344).GalacticLong = -0.446470013932275;
  vtTarget.at(345).GalacticLong = -2.2756495856466;
  vtTarget.at(346).GalacticLong = 0.837817864066572;
  vtTarget.at(347).GalacticLong = -1.33902345482271;
  vtTarget.at(348).GalacticLong = 1.17725668034362;
  vtTarget.at(349).GalacticLong = -2.84214963925724;
  vtTarget.at(350).GalacticLong = -2.31087651948803;
  vtTarget.at(351).GalacticLong = 1.37739980739096;
  vtTarget.at(352).GalacticLong = 1.70776148569054;
  vtTarget.at(353).GalacticLong = 3.14906931152967;
  vtTarget.at(354).GalacticLong = -2.60827157703617;
  vtTarget.at(355).GalacticLong = 2.8443407097;
  vtTarget.at(356).GalacticLong = -2.9074341261836;
  vtTarget.at(357).GalacticLong = 1.54874828772592;
  vtTarget.at(358).GalacticLong = -4.45223645844055;
  vtTarget.at(359).GalacticLong = 1.54311585487036;
  vtTarget.at(360).GalacticLong = -1.30030093002894;
  vtTarget.at(361).GalacticLong = -2.70238503663202;
  vtTarget.at(362).GalacticLong = 1.43118065054069;
  vtTarget.at(363).GalacticLong = 2.27891897155742;
  vtTarget.at(364).GalacticLong = -3.43501094142336;
  vtTarget.at(365).GalacticLong = 0.281110126408728;
  vtTarget.at(366).GalacticLong = -1.05230901369421;
  vtTarget.at(367).GalacticLong = -3.11626592129852;
  vtTarget.at(368).GalacticLong = -1.12621308905312;
  vtTarget.at(369).GalacticLong = 1.15587056964903;
  vtTarget.at(370).GalacticLong = -2.34604137165128;
  vtTarget.at(371).GalacticLong = -0.135582990707514;
  vtTarget.at(372).GalacticLong = 1.35909089348475;
  vtTarget.at(373).GalacticLong = -0.831217886266915;
  vtTarget.at(374).GalacticLong = -4.22335905001952;
  vtTarget.at(375).GalacticLong = 1.86059785100278;
  vtTarget.at(376).GalacticLong = 2.01966168575751;
  vtTarget.at(377).GalacticLong = 4.39247068914348;
  vtTarget.at(378).GalacticLong = -3.22061166737271;
  vtTarget.at(379).GalacticLong = -1.43273457537659;
  vtTarget.at(380).GalacticLong = -2.50992954777345;
  vtTarget.at(381).GalacticLong = 0.999397773154417;
  vtTarget.at(382).GalacticLong = 2.80476163015599;
  vtTarget.at(383).GalacticLong = -3.81592019532854;
  vtTarget.at(384).GalacticLong = -3.05669264859851;
  vtTarget.at(385).GalacticLong = -4.0073960584838;
  vtTarget.at(386).GalacticLong = 2.41840132074893;
  vtTarget.at(387).GalacticLong = -1.5148832431201;
  vtTarget.at(388).GalacticLong = 4.43443526654646;
  vtTarget.at(389).GalacticLong = -3.16800818349084;
  vtTarget.at(390).GalacticLong = 1.72905123256422;
  vtTarget.at(391).GalacticLong = -1.44952540180167;
  vtTarget.at(392).GalacticLong = 1.53143039826025;
  vtTarget.at(393).GalacticLong = 2.29772456006422;
  vtTarget.at(394).GalacticLong = -2.64235234558034;
  vtTarget.at(395).GalacticLong = -4.70264836607679;
  vtTarget.at(396).GalacticLong = 4.34380296481172;
  vtTarget.at(397).GalacticLong = -0.752432417692562;
  vtTarget.at(398).GalacticLong = 0.563067029431106;
  vtTarget.at(399).GalacticLong = 3.16913602988866;
  vtTarget.at(400).GalacticLong = 4.90854426011745;
  vtTarget.at(401).GalacticLong = 3.89467502507003;
  vtTarget.at(402).GalacticLong = 1.58030767695315;
  vtTarget.at(403).GalacticLong = 4.63640061276834;
  vtTarget.at(404).GalacticLong = -2.55997815049693;
  vtTarget.at(405).GalacticLong = 0.593819216647922;
  vtTarget.at(406).GalacticLong = -1.1056130944589;
  vtTarget.at(407).GalacticLong = 2.97221203306203;
  vtTarget.at(408).GalacticLong = -1.37917117336569;
  vtTarget.at(409).GalacticLong = 0.132470236058353;
  vtTarget.at(410).GalacticLong = -2.93586782139886;
  vtTarget.at(411).GalacticLong = 3.15995346610359;
  vtTarget.at(412).GalacticLong = -3.09954286920293;
  vtTarget.at(413).GalacticLong = 4.24642221232851;
  vtTarget.at(414).GalacticLong = 3.44857158633672;
  vtTarget.at(415).GalacticLong = 3.27328456778912;
  vtTarget.at(416).GalacticLong = -4.54901393804465;
  vtTarget.at(417).GalacticLong = -2.564523557269;
  vtTarget.at(418).GalacticLong = -2.63117955924139;
  vtTarget.at(419).GalacticLong = -3.90338255226201;
  vtTarget.at(420).GalacticLong = 4.2286044490619;
  vtTarget.at(421).GalacticLong = -4.42639666032453;
  vtTarget.at(422).GalacticLong = -1.97463365604113;
  vtTarget.at(423).GalacticLong = -0.0844412408846047;
  vtTarget.at(424).GalacticLong = -4.65414467110396;
  vtTarget.at(425).GalacticLong = -3.62458428994673;
  vtTarget.at(426).GalacticLong = -4.50044104699488;
  vtTarget.at(427).GalacticLong = -2.94266912889369;
  vtTarget.at(428).GalacticLong = -0.337114569920983;
  vtTarget.at(429).GalacticLong = -4.04672133667884;
  vtTarget.at(430).GalacticLong = -1.46223215866822;
  vtTarget.at(431).GalacticLong = -0.363462200119766;
  vtTarget.at(432).GalacticLong = -2.31890781504933;
  vtTarget.at(433).GalacticLong = -4.67700626768939;
  vtTarget.at(434).GalacticLong = 2.33822921844134;
  vtTarget.at(435).GalacticLong = 2.95865172966758;
  vtTarget.at(436).GalacticLong = 0.945039895811176;
  vtTarget.at(437).GalacticLong = -3.58081722109443;
  vtTarget.at(438).GalacticLong = 2.94769036816116;
  vtTarget.at(439).GalacticLong = -0.530877872180707;
  vtTarget.at(440).GalacticLong = -4.05398071877106;
  vtTarget.at(441).GalacticLong = -1.04518000693195;
  vtTarget.at(442).GalacticLong = 0.460226409745729;
  vtTarget.at(443).GalacticLong = 2.47292060823167;
  vtTarget.at(444).GalacticLong = 2.56744350202469;
  vtTarget.at(445).GalacticLong = 0.957628826392711;
  vtTarget.at(446).GalacticLong = 0.664171549220533;
  vtTarget.at(447).GalacticLong = -0.0818902911906694;
  vtTarget.at(448).GalacticLong = 3.15336927213071;
  vtTarget.at(449).GalacticLong = -0.340556701859795;
  vtTarget.at(450).GalacticLong = -3.41648757324352;
  vtTarget.at(451).GalacticLong = -2.96491602509415;
  vtTarget.at(452).GalacticLong = 3.12455148528506;
  vtTarget.at(453).GalacticLong = 4.00165939944644;
  vtTarget.at(454).GalacticLong = -2.83381082375235;
  vtTarget.at(455).GalacticLong = -0.845121522309077;
  vtTarget.at(456).GalacticLong = -3.87685944116607;
  vtTarget.at(457).GalacticLong = -1.02210023546687;
  vtTarget.at(458).GalacticLong = -3.4525660668889;
  vtTarget.at(459).GalacticLong = -4.39866419384027;
  vtTarget.at(460).GalacticLong = 0.327601122131169;
  vtTarget.at(461).GalacticLong = 0.248161294453325;
  vtTarget.at(462).GalacticLong = -2.45999362017802;
  vtTarget.at(463).GalacticLong = 2.74052179253864;
  vtTarget.at(464).GalacticLong = -2.51402841224668;
  vtTarget.at(465).GalacticLong = -2.44602139080542;
  vtTarget.at(466).GalacticLong = -4.54705211392832;
  vtTarget.at(467).GalacticLong = -3.04711514544885;
  vtTarget.at(468).GalacticLong = 2.16634701962929;
  vtTarget.at(469).GalacticLong = 2.69469038332509;
  vtTarget.at(470).GalacticLong = 0.961189105821534;
  vtTarget.at(471).GalacticLong = 4.48610113516692;
  vtTarget.at(472).GalacticLong = -4.43599005850884;
  vtTarget.at(473).GalacticLong = -1.34584144295821;
  vtTarget.at(474).GalacticLong = 1.48154638836362;
  vtTarget.at(475).GalacticLong = 0.0166798526503156;
  vtTarget.at(476).GalacticLong = 3.86137975219876;
  vtTarget.at(477).GalacticLong = 3.03422633653278;
  vtTarget.at(478).GalacticLong = -1.51789261917207;
  vtTarget.at(479).GalacticLong = 3.36153174752251;
  vtTarget.at(480).GalacticLong = -2.59363108129628;
  vtTarget.at(481).GalacticLong = 2.17239588052384;
  vtTarget.at(482).GalacticLong = -0.685511802856734;
  vtTarget.at(483).GalacticLong = 4.26041939496075;
  vtTarget.at(484).GalacticLong = -1.01893347430237;
  vtTarget.at(485).GalacticLong = 0.557134932875702;
  vtTarget.at(486).GalacticLong = 3.69214909323833;
  vtTarget.at(487).GalacticLong = 1.48303570121966;
  vtTarget.at(488).GalacticLong = 4.15750461809983;
  vtTarget.at(489).GalacticLong = -3.02175423019481;
  vtTarget.at(490).GalacticLong = -2.18880416570185;
  vtTarget.at(491).GalacticLong = 1.9667793115841;
  vtTarget.at(492).GalacticLong = 4.69828761728582;
  vtTarget.at(493).GalacticLong = -3.02944352870902;
  vtTarget.at(494).GalacticLong = -4.94347879142042;
  vtTarget.at(495).GalacticLong = 1.41883217695756;
  vtTarget.at(496).GalacticLong = 0.66486988334934;
  vtTarget.at(497).GalacticLong = -4.39586087206818;
  vtTarget.at(498).GalacticLong = -2.98720289635683;
  vtTarget.at(499).GalacticLong = -2.9451347974614;
  vtTarget.at(500).GalacticLong = 0.22640935854872;
  vtTarget.at(501).GalacticLong = 4.64179053273203;
  vtTarget.at(502).GalacticLong = 3.27983368432915;
  vtTarget.at(503).GalacticLong = 1.39233131357977;
  vtTarget.at(504).GalacticLong = -0.73631268373362;
  vtTarget.at(505).GalacticLong = -3.91553834856352;
  vtTarget.at(506).GalacticLong = 2.0224051213655;
  vtTarget.at(507).GalacticLong = 0.555013065294752;
  vtTarget.at(508).GalacticLong = 2.99681058650917;
  vtTarget.at(509).GalacticLong = 2.65823157002638;
  vtTarget.at(510).GalacticLong = -1.28156792108173;
  vtTarget.at(511).GalacticLong = -4.84842877059571;
  vtTarget.at(512).GalacticLong = -3.33289879508633;
  vtTarget.at(513).GalacticLong = -3.51833554754349;
  vtTarget.at(514).GalacticLong = 3.99192045189953;
  vtTarget.at(515).GalacticLong = 0.956734034237349;
  vtTarget.at(516).GalacticLong = -1.48879397140659;
  vtTarget.at(517).GalacticLong = -3.54271180684312;
  vtTarget.at(518).GalacticLong = 0.970604359320716;
  vtTarget.at(519).GalacticLong = 0.964965632129235;
  vtTarget.at(520).GalacticLong = -4.28619193550172;
  vtTarget.at(521).GalacticLong = -1.559888161427;
  vtTarget.at(522).GalacticLong = 2.95885815228912;
  vtTarget.at(523).GalacticLong = 1.61155397132166;
  vtTarget.at(524).GalacticLong = 0.841632749371045;
  vtTarget.at(525).GalacticLong = 4.39233603115224;
  vtTarget.at(526).GalacticLong = -1.49281941995577;
  vtTarget.at(527).GalacticLong = 3.96327582850898;
  vtTarget.at(528).GalacticLong = -2.3492655374272;
  vtTarget.at(529).GalacticLong = 3.31701354445718;
  vtTarget.at(530).GalacticLong = -2.23390646398202;
  vtTarget.at(531).GalacticLong = 0.406046612187303;
  vtTarget.at(532).GalacticLong = -2.03006043066617;
  vtTarget.at(533).GalacticLong = -1.59157584229803;
  vtTarget.at(534).GalacticLong = -3.28200709157339;
  vtTarget.at(535).GalacticLong = 0.141235663192605;
  vtTarget.at(536).GalacticLong = -0.379382744081046;
  vtTarget.at(537).GalacticLong = 4.77594011643098;
  vtTarget.at(538).GalacticLong = -0.0389186628253064;
  vtTarget.at(539).GalacticLong = 1.62752696984608;
  vtTarget.at(540).GalacticLong = 0.525718965616372;
  vtTarget.at(541).GalacticLong = 3.64667903437269;
  vtTarget.at(542).GalacticLong = -4.81482290733151;
  vtTarget.at(543).GalacticLong = -2.92357605924919;
  vtTarget.at(544).GalacticLong = -0.96185742791997;
  vtTarget.at(545).GalacticLong = 2.9363068537925;
  vtTarget.at(546).GalacticLong = -0.379796827440702;
  vtTarget.at(547).GalacticLong = -1.05851735978828;
  vtTarget.at(548).GalacticLong = -1.87141951165648;
  vtTarget.at(549).GalacticLong = -1.12161447569064;
  vtTarget.at(550).GalacticLong = -0.297286642335103;
  vtTarget.at(551).GalacticLong = -3.45019487600241;
  vtTarget.at(552).GalacticLong = -4.56105048925184;
  vtTarget.at(553).GalacticLong = 4.68289305878493;
  vtTarget.at(554).GalacticLong = 4.84661132483781;
  vtTarget.at(555).GalacticLong = 1.73294222948137;
  vtTarget.at(556).GalacticLong = 4.83716740839013;
  vtTarget.at(557).GalacticLong = -2.31430658969157;
  vtTarget.at(558).GalacticLong = -1.37432775522775;
  vtTarget.at(559).GalacticLong = 4.69431498309651;
  vtTarget.at(560).GalacticLong = -1.6409343640049;
  vtTarget.at(561).GalacticLong = -0.616100425217876;
  vtTarget.at(562).GalacticLong = -4.83288539091932;
  vtTarget.at(563).GalacticLong = -1.9125226053574;
  vtTarget.at(564).GalacticLong = 4.84209925332973;
  vtTarget.at(565).GalacticLong = -2.35267227879124;
  vtTarget.at(566).GalacticLong = 1.32889437596418;
  vtTarget.at(567).GalacticLong = -4.5326423722891;
  vtTarget.at(568).GalacticLong = 4.57247283145589;
  vtTarget.at(569).GalacticLong = 3.43462976527447;
  vtTarget.at(570).GalacticLong = 4.16982716417212;
  vtTarget.at(571).GalacticLong = 4.70969688255584;
  vtTarget.at(572).GalacticLong = 4.00135192300252;
  vtTarget.at(573).GalacticLong = -3.96019060450354;
  vtTarget.at(574).GalacticLong = 3.25953976662794;
  vtTarget.at(575).GalacticLong = 3.04406724970848;
  vtTarget.at(576).GalacticLong = 1.04214121784363;
  vtTarget.at(577).GalacticLong = -3.09849328978364;
  vtTarget.at(578).GalacticLong = 1.38443100743158;
  vtTarget.at(579).GalacticLong = 1.26088097178329;
  vtTarget.at(580).GalacticLong = 0.285478141915751;
  vtTarget.at(581).GalacticLong = 3.89081548292962;
  vtTarget.at(582).GalacticLong = 4.49670543861093;
  vtTarget.at(583).GalacticLong = 2.65669975414409;
  vtTarget.at(584).GalacticLong = 4.25934734693319;
  vtTarget.at(585).GalacticLong = -0.255873993496558;
  vtTarget.at(586).GalacticLong = 4.22688508754395;
  vtTarget.at(587).GalacticLong = 2.90731832079069;
  vtTarget.at(588).GalacticLong = -1.68755551574796;
  vtTarget.at(589).GalacticLong = -4.57007577809406;
  vtTarget.at(590).GalacticLong = 1.40043322942518;
  vtTarget.at(591).GalacticLong = 3.34504136549648;
  vtTarget.at(592).GalacticLong = 1.93987896243862;
  vtTarget.at(593).GalacticLong = 0.0191326436697832;
  vtTarget.at(594).GalacticLong = 1.68555770423482;
  vtTarget.at(595).GalacticLong = -1.39328309120566;
  vtTarget.at(596).GalacticLong = 1.44023414422601;
  vtTarget.at(597).GalacticLong = -0.720974097209232;
  vtTarget.at(598).GalacticLong = 1.80342722408169;
  vtTarget.at(599).GalacticLong = 2.63998586815063;
  vtTarget.at(600).GalacticLong = 0.238628172281783;
  vtTarget.at(601).GalacticLong = 0.211937301640256;
  vtTarget.at(602).GalacticLong = 1.96123472666648;
  vtTarget.at(603).GalacticLong = 3.3948317375973;
  vtTarget.at(604).GalacticLong = 0.866149229609502;
  vtTarget.at(605).GalacticLong = 0.385612779639951;
  vtTarget.at(606).GalacticLong = 1.5757999340275;
  vtTarget.at(607).GalacticLong = 3.60987843831434;
  vtTarget.at(608).GalacticLong = 2.05736464609803;
  vtTarget.at(609).GalacticLong = 1.35496859808313;
  vtTarget.at(610).GalacticLong = 1.29732251290019;
  vtTarget.at(611).GalacticLong = 2.98698402962248;
  vtTarget.at(612).GalacticLong = -3.30650198405794;
  vtTarget.at(613).GalacticLong = 2.98616322064796;
  vtTarget.at(614).GalacticLong = -3.45608982391415;
  vtTarget.at(615).GalacticLong = 1.49382900147094;
  vtTarget.at(616).GalacticLong = 0.973185276804038;
  vtTarget.at(617).GalacticLong = 3.43565503064424;
  vtTarget.at(618).GalacticLong = 1.46513469650847;
  vtTarget.at(619).GalacticLong = -3.06965612857663;
  vtTarget.at(620).GalacticLong = -2.56026111557558;
  vtTarget.at(621).GalacticLong = -1.40170557331422;
  vtTarget.at(622).GalacticLong = -3.18805286898936;
  vtTarget.at(623).GalacticLong = 2.33229074406775;
  vtTarget.at(624).GalacticLong = 4.04757672041858;
  vtTarget.at(625).GalacticLong = -0.313623071585583;
  vtTarget.at(626).GalacticLong = -3.53352418992262;
  vtTarget.at(627).GalacticLong = -1.57260992423088;
  vtTarget.at(628).GalacticLong = 0.991999875940245;
  vtTarget.at(629).GalacticLong = 1.3655912549661;
  vtTarget.at(630).GalacticLong = 2.47826636164046;
  vtTarget.at(631).GalacticLong = -0.689862489868551;
  vtTarget.at(632).GalacticLong = -2.87150652260714;
  vtTarget.at(633).GalacticLong = -1.52431198211256;
  vtTarget.at(634).GalacticLong = 3.29395689087398;
  vtTarget.at(635).GalacticLong = 3.95291688824048;
  vtTarget.at(636).GalacticLong = 3.35435835686888;
  vtTarget.at(637).GalacticLong = 4.10083468290192;
  vtTarget.at(638).GalacticLong = 0.304506606332176;
  vtTarget.at(639).GalacticLong = -1.86524473715132;
  vtTarget.at(640).GalacticLong = 3.20529328252724;
  vtTarget.at(641).GalacticLong = -2.3208848834527;
  vtTarget.at(642).GalacticLong = 1.35470183219136;
  vtTarget.at(643).GalacticLong = -3.07663280448683;
  vtTarget.at(644).GalacticLong = -1.39382877199117;
  vtTarget.at(645).GalacticLong = 0.231062724090081;
  vtTarget.at(646).GalacticLong = -0.865929749491116;
  vtTarget.at(647).GalacticLong = 2.16290948190196;
  vtTarget.at(648).GalacticLong = -2.8765109591091;
  vtTarget.at(649).GalacticLong = -2.83656386699607;
  vtTarget.at(650).GalacticLong = 1.83818328359436;
  vtTarget.at(651).GalacticLong = 2.41963578495561;
  vtTarget.at(652).GalacticLong = 4.66962055727317;
  vtTarget.at(653).GalacticLong = 2.00751593206376;
  vtTarget.at(654).GalacticLong = 1.18477289492431;
  vtTarget.at(655).GalacticLong = -1.13380023088471;
  vtTarget.at(656).GalacticLong = 0.590586670597326;
  vtTarget.at(657).GalacticLong = -4.43767028906375;
  vtTarget.at(658).GalacticLong = -1.87182900895637;
  vtTarget.at(659).GalacticLong = -0.264779309917382;
  vtTarget.at(660).GalacticLong = 0.433697547927353;
  vtTarget.at(661).GalacticLong = -4.16791403067828;
  vtTarget.at(662).GalacticLong = -0.442830864266367;
  vtTarget.at(663).GalacticLong = -4.31111221679941;
  vtTarget.at(664).GalacticLong = -4.71770037549031;
  vtTarget.at(665).GalacticLong = -4.1161941017198;
  vtTarget.at(666).GalacticLong = -3.39695508464699;
  vtTarget.at(667).GalacticLong = 0.0588430971937477;
  vtTarget.at(668).GalacticLong = 4.72783470906919;
  vtTarget.at(669).GalacticLong = 3.16780936077964;
  vtTarget.at(670).GalacticLong = -1.33254357254325;
  vtTarget.at(671).GalacticLong = -3.40422267452902;
  vtTarget.at(672).GalacticLong = 4.63643245442597;
  vtTarget.at(673).GalacticLong = 2.44234042012287;
  vtTarget.at(674).GalacticLong = -4.70129534369187;
  vtTarget.at(675).GalacticLong = -1.95410946465114;
  vtTarget.at(676).GalacticLong = -3.04809311454846;
  vtTarget.at(677).GalacticLong = -4.97136023778849;
  vtTarget.at(678).GalacticLong = -4.48363673823029;
  vtTarget.at(679).GalacticLong = -0.551394011301038;
  vtTarget.at(680).GalacticLong = 3.70847200629616;
  vtTarget.at(681).GalacticLong = 2.88503374260594;
  vtTarget.at(682).GalacticLong = -3.26604839593247;
  vtTarget.at(683).GalacticLong = 2.05130823239596;
  vtTarget.at(684).GalacticLong = -0.461495418073423;
  vtTarget.at(685).GalacticLong = 2.37379627306186;
  vtTarget.at(686).GalacticLong = 2.72250741615827;
  vtTarget.at(687).GalacticLong = -0.579837544139862;
  vtTarget.at(688).GalacticLong = 1.64750566543838;
  vtTarget.at(689).GalacticLong = -0.901045161147955;
  vtTarget.at(690).GalacticLong = 3.97470916752669;
  vtTarget.at(691).GalacticLong = 2.29843797647023;
  vtTarget.at(692).GalacticLong = -0.0957086051456759;
  vtTarget.at(693).GalacticLong = -1.22871297982984;
  vtTarget.at(694).GalacticLong = -0.371433201019408;
  vtTarget.at(695).GalacticLong = 4.4774805705538;
  vtTarget.at(696).GalacticLong = 0.0831140484702857;
  vtTarget.at(697).GalacticLong = 3.53059005287606;
  vtTarget.at(698).GalacticLong = 2.932723808289;
  vtTarget.at(699).GalacticLong = -2.68668029693889;
  vtTarget.at(700).GalacticLong = 4.78664763807094;
  vtTarget.at(701).GalacticLong = -3.51089866735432;
  vtTarget.at(702).GalacticLong = 2.41603501224275;
  vtTarget.at(703).GalacticLong = 2.72738922393048;
  vtTarget.at(704).GalacticLong = 4.70225927273596;
  vtTarget.at(705).GalacticLong = -2.63814607080604;
  vtTarget.at(706).GalacticLong = 1.21805053789736;
  vtTarget.at(707).GalacticLong = 0.673264800665403;
  vtTarget.at(708).GalacticLong = 0.543021957904271;
  vtTarget.at(709).GalacticLong = 1.2211544063266;
  vtTarget.at(710).GalacticLong = 0.245576800981224;
  vtTarget.at(711).GalacticLong = -3.34797100602286;
  vtTarget.at(712).GalacticLong = 3.00185903473324;
  vtTarget.at(713).GalacticLong = 0.914403477637393;
  vtTarget.at(714).GalacticLong = 3.55648330656667;
  vtTarget.at(715).GalacticLong = -2.44711792438496;
  vtTarget.at(716).GalacticLong = -0.809109695016764;
  vtTarget.at(717).GalacticLong = -2.9032769125773;
  vtTarget.at(718).GalacticLong = 3.6260258842331;
  vtTarget.at(719).GalacticLong = -3.11182388953935;
  vtTarget.at(720).GalacticLong = -2.25399844110657;
  vtTarget.at(721).GalacticLong = -1.09495840472136;
  vtTarget.at(722).GalacticLong = -1.07175943424775;
  vtTarget.at(723).GalacticLong = -3.94628403911998;
  vtTarget.at(724).GalacticLong = -4.87262357746592;
  vtTarget.at(725).GalacticLong = 0.898506452311258;
  vtTarget.at(726).GalacticLong = 3.16041061078526;
  vtTarget.at(727).GalacticLong = -1.83296344702443;
  vtTarget.at(728).GalacticLong = 1.48403040801313;
  vtTarget.at(729).GalacticLong = -2.42418938423296;
  vtTarget.at(730).GalacticLong = -0.619037753748151;
  vtTarget.at(731).GalacticLong = -0.867716469014244;
  vtTarget.at(732).GalacticLong = 1.11217841507172;
  vtTarget.at(733).GalacticLong = 0.606924509463834;
  vtTarget.at(734).GalacticLong = 0.805872536728623;
  vtTarget.at(735).GalacticLong = 3.63760911705304;
  vtTarget.at(736).GalacticLong = -4.24104734110036;
  vtTarget.at(737).GalacticLong = -3.17576509303902;
  vtTarget.at(738).GalacticLong = -3.95110547876183;
  vtTarget.at(739).GalacticLong = -3.02761871248911;
  vtTarget.at(740).GalacticLong = 3.16368232241493;
  vtTarget.at(741).GalacticLong = 2.77406979619376;
  vtTarget.at(742).GalacticLong = 1.41181329782659;
  vtTarget.at(743).GalacticLong = 4.39186677412417;
  vtTarget.at(744).GalacticLong = 3.96328215956836;
  vtTarget.at(745).GalacticLong = 4.0571127529779;
  vtTarget.at(746).GalacticLong = 1.62302915035213;
  vtTarget.at(747).GalacticLong = -1.66111368725968;
  vtTarget.at(748).GalacticLong = -0.862229412393924;
  vtTarget.at(749).GalacticLong = -4.58748218238412;
  vtTarget.at(750).GalacticLong = -3.30769938260372;
  vtTarget.at(751).GalacticLong = 4.59572097155832;
  vtTarget.at(752).GalacticLong = 1.79112615448143;
  vtTarget.at(753).GalacticLong = -3.69305741569545;
  vtTarget.at(754).GalacticLong = 2.55382742918972;
  vtTarget.at(755).GalacticLong = 4.04599234548232;
  vtTarget.at(756).GalacticLong = -3.67382735713279;
  vtTarget.at(757).GalacticLong = -1.06096017892518;
  vtTarget.at(758).GalacticLong = -2.4209702307184;
  vtTarget.at(759).GalacticLong = -2.07328448203102;
  vtTarget.at(760).GalacticLong = 1.81515004532967;
  vtTarget.at(761).GalacticLong = -4.55169205297958;
  vtTarget.at(762).GalacticLong = 2.8038518683555;
  vtTarget.at(763).GalacticLong = -2.00164639062803;
  vtTarget.at(764).GalacticLong = -4.63355002909262;
  vtTarget.at(765).GalacticLong = -3.01711755601858;
  vtTarget.at(766).GalacticLong = -0.646725507267814;
  vtTarget.at(767).GalacticLong = 3.26878836185829;
  vtTarget.at(768).GalacticLong = 3.06053859506795;
  vtTarget.at(769).GalacticLong = 2.39202358794731;
  vtTarget.at(770).GalacticLong = -3.67024676280738;
  vtTarget.at(771).GalacticLong = 3.87367832740206;
  vtTarget.at(772).GalacticLong = -3.79181714633684;
  vtTarget.at(773).GalacticLong = -4.25227468218517;
  vtTarget.at(774).GalacticLong = -0.00367097127898308;
  vtTarget.at(775).GalacticLong = -4.52338441356502;
  vtTarget.at(776).GalacticLong = -1.15088681080966;
  vtTarget.at(777).GalacticLong = -0.561097623399863;
  vtTarget.at(778).GalacticLong = -1.11097820139327;
  vtTarget.at(779).GalacticLong = 1.51321616618452;
  vtTarget.at(780).GalacticLong = 3.96659190299522;
  vtTarget.at(781).GalacticLong = -3.78329231573916;
  vtTarget.at(782).GalacticLong = -3.13111553410942;
  vtTarget.at(783).GalacticLong = 0.666917793621334;
  vtTarget.at(784).GalacticLong = 4.48396239719076;
  vtTarget.at(785).GalacticLong = -0.796256028530761;
  vtTarget.at(786).GalacticLong = 4.9210166459725;
  vtTarget.at(787).GalacticLong = 4.87452757389618;
  vtTarget.at(788).GalacticLong = 3.13891635877295;
  vtTarget.at(789).GalacticLong = 0.238625141029861;
  vtTarget.at(790).GalacticLong = 3.10042540871342;
  vtTarget.at(791).GalacticLong = 0.64025859069589;
  vtTarget.at(792).GalacticLong = -4.25906397813096;
  vtTarget.at(793).GalacticLong = -4.32344304727265;
  vtTarget.at(794).GalacticLong = -1.02378951240713;
  vtTarget.at(795).GalacticLong = -3.37268772852032;
  vtTarget.at(796).GalacticLong = -2.54527403635823;
  vtTarget.at(797).GalacticLong = 1.44812785637322;
  vtTarget.at(798).GalacticLong = 4.67831353190925;
  vtTarget.at(799).GalacticLong = 1.73852436133444;
  vtTarget.at(800).GalacticLong = -3.77282037662681;
  vtTarget.at(801).GalacticLong = 0.212531763708313;
  vtTarget.at(802).GalacticLong = -1.67972920394205;
  vtTarget.at(803).GalacticLong = 4.63358640874268;
  vtTarget.at(804).GalacticLong = 4.97840493444349;
  vtTarget.at(805).GalacticLong = -1.60073258127695;
  vtTarget.at(806).GalacticLong = 0.15858994644471;
  vtTarget.at(807).GalacticLong = -1.77859125237669;
  vtTarget.at(808).GalacticLong = -0.641690641671109;
  vtTarget.at(809).GalacticLong = -2.56593328145985;
  vtTarget.at(810).GalacticLong = -4.00873582719368;
  vtTarget.at(811).GalacticLong = -3.90446946923156;
  vtTarget.at(812).GalacticLong = 0.0901270632823339;
  vtTarget.at(813).GalacticLong = -2.94535784042694;
  vtTarget.at(814).GalacticLong = -0.536386837604059;
  vtTarget.at(815).GalacticLong = 1.58797303975512;
  vtTarget.at(816).GalacticLong = 3.17934914761934;
  vtTarget.at(817).GalacticLong = -2.64097679502823;
  vtTarget.at(818).GalacticLong = 0.307737602116384;
  vtTarget.at(819).GalacticLong = -4.34944816660948;
  vtTarget.at(820).GalacticLong = 3.29458620287883;
  vtTarget.at(821).GalacticLong = 4.73284242869211;
  vtTarget.at(822).GalacticLong = -3.43069138825471;
  vtTarget.at(823).GalacticLong = 0.171196594195477;
  vtTarget.at(824).GalacticLong = 3.35407768480631;
  vtTarget.at(825).GalacticLong = 0.19409010677272;
  vtTarget.at(826).GalacticLong = 4.81320262879805;
  vtTarget.at(827).GalacticLong = 1.68513926556181;
  vtTarget.at(828).GalacticLong = -4.62398652261595;
  vtTarget.at(829).GalacticLong = 4.018334023074;
  vtTarget.at(830).GalacticLong = -0.510260103411885;
  vtTarget.at(831).GalacticLong = -4.7602495914703;
  vtTarget.at(832).GalacticLong = -1.36614504006916;
  vtTarget.at(833).GalacticLong = -1.44355745126341;
  vtTarget.at(834).GalacticLong = 1.17923493444283;
  vtTarget.at(835).GalacticLong = -0.857139472141069;
  vtTarget.at(836).GalacticLong = -3.0095372230455;
  vtTarget.at(837).GalacticLong = 0.473519448835046;
  vtTarget.at(838).GalacticLong = 0.131576448787683;
  vtTarget.at(839).GalacticLong = 1.87169103822144;
  vtTarget.at(840).GalacticLong = -4.90149519818756;
  vtTarget.at(841).GalacticLong = -0.339612984784057;
  vtTarget.at(842).GalacticLong = -3.40762084232805;
  vtTarget.at(843).GalacticLong = 4.28481698572728;
  vtTarget.at(844).GalacticLong = 3.43200874687238;
  vtTarget.at(845).GalacticLong = -0.41274572787358;
  vtTarget.at(846).GalacticLong = -1.6063886187653;
  vtTarget.at(847).GalacticLong = -4.42675867882919;
  vtTarget.at(848).GalacticLong = -3.74750474658722;
  vtTarget.at(849).GalacticLong = 4.23841337855902;
  vtTarget.at(850).GalacticLong = 2.89718940295852;
  vtTarget.at(851).GalacticLong = -1.09953982969691;
  vtTarget.at(852).GalacticLong = -0.279276309268454;
  vtTarget.at(853).GalacticLong = -2.13802440808108;
  vtTarget.at(854).GalacticLong = 1.97521128572085;
  vtTarget.at(855).GalacticLong = -1.23950606224947;
  vtTarget.at(856).GalacticLong = -0.933757974816806;
  vtTarget.at(857).GalacticLong = -2.77375638438602;
  vtTarget.at(858).GalacticLong = 0.999268958596029;
  vtTarget.at(859).GalacticLong = 4.2474999766556;
  vtTarget.at(860).GalacticLong = 0.502295953144372;
  vtTarget.at(861).GalacticLong = 1.08152142347482;
  vtTarget.at(862).GalacticLong = 1.55466522250914;
  vtTarget.at(863).GalacticLong = -3.96621476643603;
  vtTarget.at(864).GalacticLong = -0.250970254859424;
  vtTarget.at(865).GalacticLong = 3.30526436659924;
  vtTarget.at(866).GalacticLong = -2.42033527672267;
  vtTarget.at(867).GalacticLong = 4.41942612840803;
  vtTarget.at(868).GalacticLong = -1.70686062050509;
  vtTarget.at(869).GalacticLong = -3.81490801212563;
  vtTarget.at(870).GalacticLong = 1.02715160787824;
  vtTarget.at(871).GalacticLong = 0.110825668790886;
  vtTarget.at(872).GalacticLong = 1.70129432033198;
  vtTarget.at(873).GalacticLong = 4.85252418577305;
  vtTarget.at(874).GalacticLong = 2.97696466647887;
  vtTarget.at(875).GalacticLong = 0.635425848513373;
  vtTarget.at(876).GalacticLong = 3.34925795389925;
  vtTarget.at(877).GalacticLong = 3.53307078382121;
  vtTarget.at(878).GalacticLong = -1.52143362506974;
  vtTarget.at(879).GalacticLong = -2.94425317954483;
  vtTarget.at(880).GalacticLong = -3.35053018344137;
  vtTarget.at(881).GalacticLong = 4.72641617861745;
  vtTarget.at(882).GalacticLong = 2.5222226555056;
  vtTarget.at(883).GalacticLong = 4.83635627990889;
  vtTarget.at(884).GalacticLong = 3.03592811667383;
  vtTarget.at(885).GalacticLong = 2.48885794673877;
  vtTarget.at(886).GalacticLong = 3.35925915645898;
  vtTarget.at(887).GalacticLong = 4.95675515125096;
  vtTarget.at(888).GalacticLong = 0.0413947413390723;
  vtTarget.at(889).GalacticLong = -3.81838940878467;
  vtTarget.at(890).GalacticLong = 4.90577822435938;
  vtTarget.at(891).GalacticLong = -2.67756021708913;
  vtTarget.at(892).GalacticLong = 1.59349513358035;
  vtTarget.at(893).GalacticLong = -4.44078658484419;
  vtTarget.at(894).GalacticLong = -3.33808467714078;
  vtTarget.at(895).GalacticLong = -3.44117671898985;
  vtTarget.at(896).GalacticLong = -4.45663075299566;
  vtTarget.at(897).GalacticLong = 0.776684583948386;
  vtTarget.at(898).GalacticLong = -0.199801092487769;
  vtTarget.at(899).GalacticLong = 4.98337052950952;
  vtTarget.at(900).GalacticLong = -0.279500326935787;
  vtTarget.at(901).GalacticLong = 3.42407406968367;
  vtTarget.at(902).GalacticLong = 2.9811780883429;
  vtTarget.at(903).GalacticLong = -3.21804771598297;
  vtTarget.at(904).GalacticLong = -3.57195037559013;
  vtTarget.at(905).GalacticLong = -4.46241634667767;
  vtTarget.at(906).GalacticLong = 2.91103610134786;
  vtTarget.at(907).GalacticLong = -1.86486720054109;
  vtTarget.at(908).GalacticLong = -3.32788664187672;
  vtTarget.at(909).GalacticLong = -0.0128261125181268;
  vtTarget.at(910).GalacticLong = -0.587718478251413;
  vtTarget.at(911).GalacticLong = -0.97669647620926;
  vtTarget.at(912).GalacticLong = -2.08679767574076;
  vtTarget.at(913).GalacticLong = 1.96963541608032;
  vtTarget.at(914).GalacticLong = 1.42774305353651;
  vtTarget.at(915).GalacticLong = 0.456147397654377;
  vtTarget.at(916).GalacticLong = -1.72904360122941;
  vtTarget.at(917).GalacticLong = 1.35951396990348;
  vtTarget.at(918).GalacticLong = -1.70861491352497;
  vtTarget.at(919).GalacticLong = -2.39346162030893;
  vtTarget.at(920).GalacticLong = 1.84787791484694;
  vtTarget.at(921).GalacticLong = -0.305546109032115;
  vtTarget.at(922).GalacticLong = 1.99053079534548;
  vtTarget.at(923).GalacticLong = 2.52294906104718;
  vtTarget.at(924).GalacticLong = 1.8801055497586;
  vtTarget.at(925).GalacticLong = 0.187565617856235;
  vtTarget.at(926).GalacticLong = -3.43768994394042;
  vtTarget.at(927).GalacticLong = -1.8873948379779;
  vtTarget.at(928).GalacticLong = -4.6948056029357;
  vtTarget.at(929).GalacticLong = 2.55968727980534;
  vtTarget.at(930).GalacticLong = 0.196877649734111;
  vtTarget.at(931).GalacticLong = -1.1154200849335;
  vtTarget.at(932).GalacticLong = -1.64574341589614;
  vtTarget.at(933).GalacticLong = 4.59212091989645;
  vtTarget.at(934).GalacticLong = -0.838199518415067;
  vtTarget.at(935).GalacticLong = -3.5863894491427;
  vtTarget.at(936).GalacticLong = 0.655197226748675;
  vtTarget.at(937).GalacticLong = 1.92667383750826;
  vtTarget.at(938).GalacticLong = 1.01952714656663;
  vtTarget.at(939).GalacticLong = -2.35006217670683;
  vtTarget.at(940).GalacticLong = -3.07224929232082;
  vtTarget.at(941).GalacticLong = -4.93084655761927;
  vtTarget.at(942).GalacticLong = 0.73075963452739;
  vtTarget.at(943).GalacticLong = 1.54332771695762;
  vtTarget.at(944).GalacticLong = -4.75169718659088;
  vtTarget.at(945).GalacticLong = 2.92859800930729;
  vtTarget.at(946).GalacticLong = 4.92272432457932;
  vtTarget.at(947).GalacticLong = 3.52373296910141;
  vtTarget.at(948).GalacticLong = 3.83136975513142;
  vtTarget.at(949).GalacticLong = -1.89064682108297;
  vtTarget.at(950).GalacticLong = -3.68678577046509;
  vtTarget.at(951).GalacticLong = -4.07536792675878;
  vtTarget.at(952).GalacticLong = -4.75100364485346;
  vtTarget.at(953).GalacticLong = 0.769191410746162;
  vtTarget.at(954).GalacticLong = -3.59621459994052;
  vtTarget.at(955).GalacticLong = 3.0605502730693;
  vtTarget.at(956).GalacticLong = -1.39476991868501;
  vtTarget.at(957).GalacticLong = -4.76563805661013;
  vtTarget.at(958).GalacticLong = 0.578277444612918;
  vtTarget.at(959).GalacticLong = 3.30427042243328;
  vtTarget.at(960).GalacticLong = -0.522033985401239;
  vtTarget.at(961).GalacticLong = 0.710049553718896;
  vtTarget.at(962).GalacticLong = 4.23678320245367;
  vtTarget.at(963).GalacticLong = -0.229927341706134;
  vtTarget.at(964).GalacticLong = -1.72841286949634;
  vtTarget.at(965).GalacticLong = -3.7472463519774;
  vtTarget.at(966).GalacticLong = -1.83779208762961;
  vtTarget.at(967).GalacticLong = 2.68068667660445;
  vtTarget.at(968).GalacticLong = -2.94687365056025;
  vtTarget.at(969).GalacticLong = 1.66627131015281;
  vtTarget.at(970).GalacticLong = -0.843242592546955;
  vtTarget.at(971).GalacticLong = 2.7639088184175;
  vtTarget.at(972).GalacticLong = 0.447171951429705;
  vtTarget.at(973).GalacticLong = -1.53595550042545;
  vtTarget.at(974).GalacticLong = 0.893140530375742;
  vtTarget.at(975).GalacticLong = -1.90457997422892;
  vtTarget.at(976).GalacticLong = 3.62254174955805;
  vtTarget.at(977).GalacticLong = 1.96127167817055;
  vtTarget.at(978).GalacticLong = 3.28960563475235;
  vtTarget.at(979).GalacticLong = 1.46988238680454;
  vtTarget.at(980).GalacticLong = -3.60025044736825;
  vtTarget.at(981).GalacticLong = -3.3032685364189;
  vtTarget.at(982).GalacticLong = -2.05031566670963;
  vtTarget.at(983).GalacticLong = 1.38965731393779;
  vtTarget.at(984).GalacticLong = 3.95306543552687;
  vtTarget.at(985).GalacticLong = -2.10438716055755;
  vtTarget.at(986).GalacticLong = 0.160171298795793;
  vtTarget.at(987).GalacticLong = 1.16815637976279;
  vtTarget.at(988).GalacticLong = 1.04654975657193;
  vtTarget.at(989).GalacticLong = -2.1959563756659;
  vtTarget.at(990).GalacticLong = -2.72356638667002;
  vtTarget.at(991).GalacticLong = -4.11335690158507;
  vtTarget.at(992).GalacticLong = -0.701659036555658;
  vtTarget.at(993).GalacticLong = 1.31725322328562;
  vtTarget.at(994).GalacticLong = -2.79953251525685;
  vtTarget.at(995).GalacticLong = -0.859953667965695;
  vtTarget.at(996).GalacticLong = 0.551107528388364;
  vtTarget.at(997).GalacticLong = -2.70793897349847;
  vtTarget.at(998).GalacticLong = 2.48416538326394;
  vtTarget.at(999).GalacticLong = 2.49809889345878;
  vtTarget.at(1000).GalacticLong = 3.10554237987262;
  vtTarget.at(1001).GalacticLong = -1.5951371466996;
  vtTarget.at(1002).GalacticLong = -0.659780016972269;
  vtTarget.at(1003).GalacticLong = -0.144110572212732;
  vtTarget.at(1004).GalacticLong = 2.56583759852638;
  vtTarget.at(1005).GalacticLong = -3.53308857775611;
  vtTarget.at(1006).GalacticLong = 2.58142038849156;
  vtTarget.at(1007).GalacticLong = -3.93072768818573;
  vtTarget.at(1008).GalacticLong = 3.69769230973197;
  vtTarget.at(1009).GalacticLong = -2.52451723428768;
  vtTarget.at(1010).GalacticLong = 0.398785406835334;
  vtTarget.at(1011).GalacticLong = -2.75556419775869;
  vtTarget.at(1012).GalacticLong = 1.39737708868656;
  vtTarget.at(1013).GalacticLong = 0.164316378440357;
  vtTarget.at(1014).GalacticLong = -0.222267307278894;
  vtTarget.at(1015).GalacticLong = -0.192101581261701;
  vtTarget.at(1016).GalacticLong = -4.02538494782811;
  vtTarget.at(1017).GalacticLong = -1.13541153546498;
  vtTarget.at(1018).GalacticLong = 1.21772733628181;
  vtTarget.at(1019).GalacticLong = 4.93640885853071;
  vtTarget.at(1020).GalacticLong = -1.36292990537467;
  vtTarget.at(1021).GalacticLong = 4.43736542680526;
  vtTarget.at(1022).GalacticLong = -4.96061548286042;
  vtTarget.at(1023).GalacticLong = -3.5988125645647;
  vtTarget.at(1024).GalacticLong = -4.91810130536813;
  vtTarget.at(1025).GalacticLong = -3.70057865367723;
  vtTarget.at(1026).GalacticLong = -2.5192115982069;
  vtTarget.at(1027).GalacticLong = 3.22157346416908;
  vtTarget.at(1028).GalacticLong = -0.64898220830821;
  vtTarget.at(1029).GalacticLong = -3.9570453859326;
  vtTarget.at(1030).GalacticLong = 4.15970190512497;
  vtTarget.at(1031).GalacticLong = -2.45501150600342;
  vtTarget.at(1032).GalacticLong = 2.45939541197193;
  vtTarget.at(1033).GalacticLong = 4.57051302048394;
  vtTarget.at(1034).GalacticLong = -4.29859227523893;
  vtTarget.at(1035).GalacticLong = -1.76292316603277;
  vtTarget.at(1036).GalacticLong = -2.99013181754888;
  vtTarget.at(1037).GalacticLong = 4.08118286433693;
  vtTarget.at(1038).GalacticLong = 2.90352874853436;
  vtTarget.at(1039).GalacticLong = 0.437264452896121;
  vtTarget.at(1040).GalacticLong = -0.0522732793506329;
  vtTarget.at(1041).GalacticLong = 2.13563678019069;
  vtTarget.at(1042).GalacticLong = 2.00769806262436;
  vtTarget.at(1043).GalacticLong = 4.0759225850828;
  vtTarget.at(1044).GalacticLong = -4.84530192689921;
  vtTarget.at(1045).GalacticLong = -1.46297701708644;
  vtTarget.at(1046).GalacticLong = 3.01058077022056;
  vtTarget.at(1047).GalacticLong = -2.71930783822015;
  vtTarget.at(1048).GalacticLong = 3.69310747936147;
  vtTarget.at(1049).GalacticLong = 2.3307197796918;
  vtTarget.at(1050).GalacticLong = 3.35118290045832;
  vtTarget.at(1051).GalacticLong = 2.77071815223651;
  vtTarget.at(1052).GalacticLong = -3.65408496459147;
  vtTarget.at(1053).GalacticLong = 1.91895347039463;
  vtTarget.at(1054).GalacticLong = -0.317264876512847;
  vtTarget.at(1055).GalacticLong = -3.54559762146836;
  vtTarget.at(1056).GalacticLong = -2.76854642720173;
  vtTarget.at(1057).GalacticLong = -4.51322789949601;
  vtTarget.at(1058).GalacticLong = -1.92119257392519;
  vtTarget.at(1059).GalacticLong = -4.42485244247821;
  vtTarget.at(1060).GalacticLong = -1.08245816398576;
  vtTarget.at(1061).GalacticLong = -1.94246244056693;
  vtTarget.at(1062).GalacticLong = 0.0935360203597995;
  vtTarget.at(1063).GalacticLong = 2.12355667618791;
  vtTarget.at(1064).GalacticLong = 1.09589652031875;
  vtTarget.at(1065).GalacticLong = -2.1752693095938;
  vtTarget.at(1066).GalacticLong = 3.11360353801616;
  vtTarget.at(1067).GalacticLong = 0.573782671584965;
  vtTarget.at(1068).GalacticLong = -0.924259108982945;
  vtTarget.at(1069).GalacticLong = -3.93382501707181;
  vtTarget.at(1070).GalacticLong = -1.13217353600193;
  vtTarget.at(1071).GalacticLong = -3.22324840781758;
  vtTarget.at(1072).GalacticLong = -3.69927117803764;
  vtTarget.at(1073).GalacticLong = 3.76211730087434;
  vtTarget.at(1074).GalacticLong = -0.302314089523454;
  vtTarget.at(1075).GalacticLong = 1.73365036441965;
  vtTarget.at(1076).GalacticLong = 0.619290791417586;
  vtTarget.at(1077).GalacticLong = 4.19985668317688;
  vtTarget.at(1078).GalacticLong = 3.10047442560291;
  vtTarget.at(1079).GalacticLong = 3.67055771154627;
  vtTarget.at(1080).GalacticLong = 3.26182795300976;
  vtTarget.at(1081).GalacticLong = -4.12317123654376;
  vtTarget.at(1082).GalacticLong = 0.724876552631746;
  vtTarget.at(1083).GalacticLong = -2.30718264484;
  vtTarget.at(1084).GalacticLong = -1.22968110506546;
  vtTarget.at(1085).GalacticLong = 1.50558276680824;
  vtTarget.at(1086).GalacticLong = 3.72190887299679;
  vtTarget.at(1087).GalacticLong = -2.80714876310664;
  vtTarget.at(1088).GalacticLong = -1.66556242303171;
  vtTarget.at(1089).GalacticLong = -4.75798850626151;
  vtTarget.at(1090).GalacticLong = -3.94061740603494;
  vtTarget.at(1091).GalacticLong = -4.20989577844988;
  vtTarget.at(1092).GalacticLong = -2.04856909413064;
  vtTarget.at(1093).GalacticLong = 4.65581190561721;
  vtTarget.at(1094).GalacticLong = -1.74036447368682;
  vtTarget.at(1095).GalacticLong = -1.41210934186784;
  vtTarget.at(1096).GalacticLong = -0.252386050188576;
  vtTarget.at(1097).GalacticLong = -1.98234150639792;
  vtTarget.at(1098).GalacticLong = 3.80657706873784;
  vtTarget.at(1099).GalacticLong = -3.30219990496273;

  // Initialise variables.
  CoveringCircleFiller.GalacticLat = 0.0;
  CoveringCircleFiller.GalacticLong = 0.0;
  CoveringCircleFiller.TargetsCovered = -1;

  // Initialise MPI.
  MPI_Init(&argc, &argv);

  // Get the number of processes available in the communicator.
  MPI_Comm_size(MPI_COMM_WORLD, &intSize);
  // Get a process's rank within the communicator, from 0 to (intSize - 1).
  MPI_Comm_rank(MPI_COMM_WORLD, &intRank);

  dblTimeStarted = MPI_Wtime();

  // Create an MPI struct for the CoveringCircle struct.
  intItems = std::size(intBlocklengthsArray);
  offsetsArray[0] = offsetof(CoveringCircle, GalacticLat);
  offsetsArray[1] = offsetof(CoveringCircle, GalacticLong);
  offsetsArray[2] = offsetof(CoveringCircle, TargetsCovered);
  offsetsArray[3] = offsetof(CoveringCircle, StableID);
  MPI_Type_create_struct(intItems, intBlocklengthsArray, offsetsArray,
                         typesArray, &MPI_CoveringCirclePrototype);

  // Resize the MPI struct.
  // Ensure its length matches its length implemented in memory.
  MPI_Type_get_extent(MPI_CoveringCirclePrototype, &offsetLowerBound,
                      &offsetExtent);

  // Calculate the vector's inter-element distance.
  MPI_Get_address(&CoveringCircleExtentProbe[0], &a0);
  MPI_Get_address(&CoveringCircleExtentProbe[1], &a1);
  offsetExtent = MPI_Aint_diff(a1, a0);

  MPI_Type_create_resized(MPI_CoveringCirclePrototype, offsetLowerBound,
                          offsetExtent, &MPI_CoveringCircle);
  MPI_Type_commit(&MPI_CoveringCircle);

  do {
    // Build target pair indices array.
    for (intI = 0; intI <= (intNumTargets - 2); intI++) {
      for (intJ = (intI + 1); intJ <= (intNumTargets - 1); intJ++) {
        vtIndexPair.push_back(
            {0.0, 0.0, static_cast<int16_t>(intI), static_cast<int16_t>(intJ)});
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
    vtCoveringCircle.clear();
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

      if (isEpsilonClose(dblSeparation, FOV_DIAMETER)) {
        // Add one cover for two targets.
        vtCoveringCircle.resize(i32Covers + 1);

        // Calculate the midpoint of the chord between the targets.
        vtCoveringCircle.at(i32Covers).GalacticLat =
            0.5 * (vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLat +
                   vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLat);
        vtCoveringCircle.at(i32Covers).GalacticLong =
            0.5 * (vtTarget.at(vtIndexPair.at(intX).TargetOne).GalacticLong +
                   vtTarget.at(vtIndexPair.at(intX).TargetTwo).GalacticLong);
        // Update the cover ID.
        // N.B. "<<" is a bitwise left shift.
        vtCoveringCircle.at(i32Covers).StableID =
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetOne) << 17) |
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetTwo) << 1);

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
        vtCoveringCircle.resize(i32Covers + 1);
        // Update latitude.
        vtCoveringCircle.at(i32Covers).GalacticLat =
            dblAvgLat + dblScaledUnitPerpVectorLong;
        // Update longitude.
        vtCoveringCircle.at(i32Covers).GalacticLong =
            dblAvgLong - dblScaledUnitPerpVectorLat;
        // Update the cover ID.
        ui64StableID1 =
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetOne) << 17) |
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetTwo) << 1) | 0;
        vtCoveringCircle.at(i32Covers).StableID = ui64StableID1;
        // Increase cover counter.
        i32Covers += 1;

        // Cover 2.
        // N.B. Keep the comments in this section, or autoformatters mangle it.
        vtCoveringCircle.resize(i32Covers + 1);
        // Update latitude.
        vtCoveringCircle.at(i32Covers).GalacticLat =
            dblAvgLat - dblScaledUnitPerpVectorLong;
        // Update longitude.
        vtCoveringCircle.at(i32Covers).GalacticLong =
            dblAvgLong + dblScaledUnitPerpVectorLat;
        // Update the cover ID.
        ui64StableID12 =
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetOne) << 17) |
            (static_cast<uint64_t>(vtIndexPair.at(intX).TargetTwo) << 1) | 1;
        vtCoveringCircle.at(i32Covers).StableID = ui64StableID12;
        // Increase cover counter.
        i32Covers += 1;
      }
    }

    // Cover singleton targets.
    for (intX = 0; intX < intNumTargets; intX++) {
      vtCoveringCircle.resize(i32Covers + 1);
      // Update latitude.
      vtCoveringCircle.at(i32Covers).GalacticLat =
          vtTarget.at(intX).GalacticLat;
      // Update longitude.
      vtCoveringCircle.at(i32Covers).GalacticLong =
          vtTarget.at(intX).GalacticLong;
      // Update the cover ID.
      vtCoveringCircle.at(i32Covers).StableID =
          (1ULL << 63) | static_cast<uint64_t>(intX);
      // Increase cover counter.
      i32Covers += 1;
    }

    // Free unused local index pair data.
    vtIndexPair.clear();

    // Determine the local maximum number of covering circles created.
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

    // Give each covering circle a uniform extent to facilitate Allgather.
    if (static_cast<std::size_t>(vtCoveringCircle.size()) < i32MaxCovers) {
      vtCoveringCircle.resize(i32MaxCovers, CoveringCircleFiller);
    }

    // Build vector to receive global covering circle data into.
    vtGlobalCoveringCircle.resize(i32MaxCovers * intSize);

    // Distribute local covering circle data across all ranks.
    // N.B. ".data()" returns a pointer to the underlying raw array used by the
    // vectors to store their data.
    // N.B. MPI_Allgather is a synchronising collective communication function.
    // It doesn't need leading or trailing barriers.
    MPI_Allgather(vtCoveringCircle.data(), i32MaxCovers, MPI_CoveringCircle,
                  vtGlobalCoveringCircle.data(), i32MaxCovers,
                  MPI_CoveringCircle, MPI_COMM_WORLD);

    // Free unused covering circle data.
    vtCoveringCircle.clear();

    // Remove filler covering circles.
    for (intI = (vtGlobalCoveringCircle.size() - 1); intI >= 0; intI--) {
      if (vtGlobalCoveringCircle.at(intI).TargetsCovered == -1)
        vtGlobalCoveringCircle.erase(vtGlobalCoveringCircle.begin() + intI);
    }
    i32TotalCovers = vtGlobalCoveringCircle.size();

    // Configure parallelisation of target coverage calculations.
    GetWorkRange(i32TotalCovers, intRank, intSize, intStart, intEnd);

    // Prepare for deterministic greedy set cover tie-breaking.
    // A "stable ID" is used for consistent cover selection, regardless of
    // rank count or race conditions.
    ui64LocalBestStableID = UINT64_MAX;
    intMaxCoverIndex = -1;

    // Reset cover counts.
    for (intI = 0; intI < i32TotalCovers; intI++) {
      vtGlobalCoveringCircle.at(intI).TargetsCovered = 0;
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
                           vtGlobalCoveringCircle.at(intI).GalacticLat;
        dblSeparationLong = vtTarget.at(intJ).GalacticLong -
                            vtGlobalCoveringCircle.at(intI).GalacticLong;

        // A target is covered if it's within FOV_RADIUS of the cover's centre.
        // N.B. Use squared-distance comparison for performance:
        // Slower: if (sqrt(a*a + b*b) < FOV_RADIUS)
        // Faster: if ((a*a + b*b) < FOV_RADIUS_POW2)
        if ((dblSeparationLat * dblSeparationLat +
             dblSeparationLong * dblSeparationLong) < FOV_RADIUS_POW2) {
          vtGlobalCoveringCircle.at(intI).TargetsCovered += 1;
          intTargetsTemp += 1;
        }
      }
      // Assess greediest covers, with robust tie-breaking.
      ui64CurrentID = vtGlobalCoveringCircle.at(intI).StableID;
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

    for (const auto &GreedyCoverCandidate : vtGreedyCovers) {
      if (GreedyCoverCandidate.TargetsCovered >
              GreedyCoverGlobal.TargetsCovered ||
          (GreedyCoverCandidate.TargetsCovered ==
               GreedyCoverGlobal.TargetsCovered &&
           GreedyCoverCandidate.StableID < GreedyCoverGlobal.StableID)) {
        GreedyCoverGlobal = GreedyCoverCandidate;
      }
    }

    ui64GlobalStableID = GreedyCoverGlobal.StableID;
    intGlobalMaxCoverIndex = -1;
    for (intI = 0; intI < i32TotalCovers; intI++) {
      if (vtGlobalCoveringCircle.at(intI).StableID == ui64GlobalStableID) {
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
          vtGlobalCoveringCircle.at(intGlobalMaxCoverIndex).GalacticLat;
      dblSeparationLong =
          vtTarget.at(intJ).GalacticLong -
          vtGlobalCoveringCircle.at(intGlobalMaxCoverIndex).GalacticLong;

      // A target is covered if it's within FOV_RADIUS of the cover's centre.
      // N.B. Use squared-distance comparison for performance:
      // Slower: if (sqrt(a*a + b*b) < FOV_RADIUS)
      // Faster: if ((a*a + b*b) < FOV_RADIUS_POW2)
      if ((dblSeparationLat * dblSeparationLat +
           dblSeparationLong * dblSeparationLong) < FOV_RADIUS_POW2) {
        vtTarget.at(intJ).Covered = true;
      }
    }

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
           << vtGlobalCoveringCircle.at(intGlobalMaxCoverIndex).GalacticLat
           << endl;
      cout << "   Galactic longitude: "
           << vtGlobalCoveringCircle.at(intGlobalMaxCoverIndex).GalacticLong
           << endl;
      cout << "   Targets covered: " << intTargetsMax << endl;
    }

    // Free vector data.
    vtGlobalCoveringCircle.clear();
  } while (intNumTargets > 0);

  MPI_Barrier(MPI_COMM_WORLD);
  dblTimeCompleted = MPI_Wtime();

  // Terminate the MPI execution environment.
  MPI_Finalize();

  cout << "MPI Execution time: " << dblTimeCompleted - dblTimeStarted << endl;

  return 0;
}

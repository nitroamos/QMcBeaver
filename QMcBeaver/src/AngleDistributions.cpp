#include "AngleDistributions.h"

#define PI 3.14159265359

Array1D<double> AngleDistributions::getPhiArray(int index)
{
  Array1D<double> phi_array(21);

  if (index == 3)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.327655967106;
      phi_array(2) = 0.560578576795;
      phi_array(3) = 0.750594366479;
      phi_array(4) = 0.916377761261;
      phi_array(5) = 1.0670557307;
      phi_array(6) = 1.20795532873;
      phi_array(7) = 1.3426062987;
      phi_array(8) = 1.47362013795;
      phi_array(9) = 1.60314770425;
      phi_array(10) = 1.7331658919;
      phi_array(11) = 1.86570815646;
      phi_array(12) = 2.00311417765;
      phi_array(13) = 2.14837732602;
      phi_array(14) = 2.30577228193;
      phi_array(15) = 2.4821488742;
      phi_array(16) = 2.69021721276;
      phi_array(17) = 2.95946679675;
      phi_array(18) = 3.39907010765;
      phi_array(19) = 5.50481419047;
      phi_array(20) = 2*PI;
    }

  else if (index == 4)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.778371116713;
      phi_array(2) = 2.88411519953;
      phi_array(3) = 3.32371851043;
      phi_array(4) = 3.59296809442;
      phi_array(5) = 3.80103643298;
      phi_array(6) = 3.97741302524;
      phi_array(7) = 4.13480798116;
      phi_array(8) = 4.28007112953;
      phi_array(9) = 4.41747715072;
      phi_array(10) = 4.55001941528;
      phi_array(11) = 4.68003760293;
      phi_array(12) = 4.80956516923;
      phi_array(13) = 4.94057900848;
      phi_array(14) = 5.07522997845;
      phi_array(15) = 5.21612957648;
      phi_array(16) = 5.36680754592;
      phi_array(17) = 5.5325909407;
      phi_array(18) = 5.72260673038;
      phi_array(19) = 5.95552934008;
      phi_array(20) = 2*PI;
    } 

  else if (index == 6)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.128341746125;
      phi_array(2) = 0.258754953877;
      phi_array(3) = 0.393571324756;
      phi_array(4) = 0.535727491455;
      phi_array(5) = 0.689369470051;
      phi_array(6) = 0.861083099716;
      phi_array(7) = 1.06300466921;
      phi_array(8) = 1.32291249526;
      phi_array(9) = 1.73889905857;
      phi_array(10) = 3.14159265359;
      phi_array(11) = 4.54428624861;
      phi_array(12) = 4.96027281192;
      phi_array(13) = 5.22018063797;
      phi_array(14) = 5.42210220747;
      phi_array(15) = 5.59381583713;
      phi_array(16) = 5.74745781572;
      phi_array(17) = 5.88961398242;
      phi_array(18) = 6.02443035331;
      phi_array(19) = 6.15484356106;
      phi_array(20) = 2*PI;
    }     

  else if (index == 7)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 1.40269359503;
      phi_array(2) = 1.81868015833;
      phi_array(3) = 2.07858798438;
      phi_array(4) = 2.28050955388;
      phi_array(5) = 2.45222318354;
      phi_array(6) = 2.60586516213;
      phi_array(7) = 2.74802132883;
      phi_array(8) = 2.88283769972;
      phi_array(9) = 3.01325090747;
      phi_array(10) = 3.14159265359;
      phi_array(11) = 3.26993439971;
      phi_array(12) = 3.40034760746;
      phi_array(13) = 3.53516397834;
      phi_array(14) = 3.67732014505;
      phi_array(15) = 3.83096212364;
      phi_array(16) = 4.0026757533;
      phi_array(17) = 4.2045973228;
      phi_array(18) = 4.46450514885;
      phi_array(19) = 4.88049171216;
      phi_array(20) = 2*PI;
    }

  else if (index == 8)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.360502118928;
      phi_array(2) = 0.588676739245;
      phi_array(3) = 0.768113525124;
      phi_array(4) = 0.921862717029;
      phi_array(5) = 1.06004830735;
      phi_array(6) = 1.18820902164;
      phi_array(7) = 1.30983978629;
      phi_array(8) = 1.42742112219;
      phi_array(9) = 1.54289928979;
      phi_array(10) = 1.6579779379;
      phi_array(11) = 1.77431595382;
      phi_array(12) = 1.89371362856;
      phi_array(13) = 2.01834673471;
      phi_array(14) = 2.15113607751;
      phi_array(15) = 2.29644342412;
      phi_array(16) = 2.46166359701;
      phi_array(17) = 2.66171947816;
      phi_array(18) = 2.93741072021;
      phi_array(19) = 3.59784745921;
      phi_array(20) = 2*PI;
    }

  else if (index == 9)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.139403980084;
      phi_array(2) = 0.295082016953;
      phi_array(3) = 0.477876057221;
      phi_array(4) = 0.713191635944;
      phi_array(5) = 1.1009648807;
      phi_array(6) = 3.87475025655;
      phi_array(7) = 4.4256348224;
      phi_array(8) = 4.69026854925;
      phi_array(9) = 4.8858783257;
      phi_array(10) = 5.04870358345;
      phi_array(11) = 5.19254522274;
      phi_array(12) = 5.3243860571;
      phi_array(13) = 5.44840256003;
      phi_array(14) = 5.56741743594;
      phi_array(15) = 5.68355797714;
      phi_array(16) = 5.79860364476;
      phi_array(17) = 5.91420873452;
      phi_array(18) = 6.03208935068;
      phi_array(19) = 6.15423396889;
      phi_array(20) = 2*PI;
    }

  else if (index == 10)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 2.16510311178;
      phi_array(2) = 2.49299538786;
      phi_array(3) = 2.71150073558;
      phi_array(4) = 2.88623092642;
      phi_array(5) = 3.03723306725;
      phi_array(6) = 3.17367728436;
      phi_array(7) = 3.30069688015;
      phi_array(8) = 3.42159958177;
      phi_array(9) = 3.53876085884;
      phi_array(10) = 3.65407772619;
      phi_array(11) = 3.7692354223;
      phi_array(12) = 3.88590196091;
      phi_array(13) = 4.00591560999;
      phi_array(14) = 4.13153348047;
      phi_array(15) = 4.26582717861;
      phi_array(16) = 4.41347624058;
      phi_array(17) = 4.58256475326;
      phi_array(18) = 4.7899613285;
      phi_array(19) = 5.08509506874;
      phi_array(20) = 2*PI;
    }

  else if (index == 11)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.165243831658;
      phi_array(2) = 0.314701677956;
      phi_array(3) = 0.454153021336;
      phi_array(4) = 0.587376785657;
      phi_array(5) = 0.71715235498;
      phi_array(6) = 0.845782641838;
      phi_array(7) = 0.975424813317;
      phi_array(8) = 1.10836574495;
      phi_array(9) = 1.24734370633;
      phi_array(10) = 1.39605073367;
      phi_array(11) = 1.56010666893;
      phi_array(12) = 1.7493159641;
      phi_array(13) = 1.98436091633;
      phi_array(14) = 2.32585009432;
      phi_array(15) = 3.14159265359;
      phi_array(16) = 4.62117648083;
      phi_array(17) = 5.5001392095;
      phi_array(18) = 5.85288409818;
      phi_array(19) = 6.0919436807;
      phi_array(20) = 2*PI;
    }

  else if (index == 12)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 1.47958382724;
      phi_array(2) = 2.35854655591;
      phi_array(3) = 2.71129144459;
      phi_array(4) = 2.95035102711;
      phi_array(5) = 3.14159265359;
      phi_array(6) = 3.30683648525;
      phi_array(7) = 3.45629433155;
      phi_array(8) = 3.59574567493;
      phi_array(9) = 3.72896943925;
      phi_array(10) = 3.85874500857;
      phi_array(11) = 3.98737529543;
      phi_array(12) = 4.11701746691;
      phi_array(13) = 4.24995839854;
      phi_array(14) = 4.38893635992;
      phi_array(15) = 4.53764338726;
      phi_array(16) = 4.70169932252;
      phi_array(17) = 4.89090861769;
      phi_array(18) = 5.12595356992;
      phi_array(19) = 5.46744274791;
      phi_array(20) = 2*PI;
    }

  else if (index == 13)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.191241626477;
      phi_array(2) = 0.430301209002;
      phi_array(3) = 0.783046097674;
      phi_array(4) = 1.66200882635;
      phi_array(5) = 3.14159265359;
      phi_array(6) = 3.95733521286;
      phi_array(7) = 4.29882439085;
      phi_array(8) = 4.53386934308;
      phi_array(9) = 4.72307863825;
      phi_array(10) = 4.88713457351;
      phi_array(11) = 5.03584160085;
      phi_array(12) = 5.17481956223;
      phi_array(13) = 5.30776049386;
      phi_array(14) = 5.43740266534;
      phi_array(15) = 5.5660329522;
      phi_array(16) = 5.69580852152;
      phi_array(17) = 5.82903228584;
      phi_array(18) = 5.96848362922;
      phi_array(19) = 6.11794147552;
      phi_array(20) = 2*PI;
    }
      
  else if (index == 14)
    {
      phi_array(0) = 0.0;
      phi_array(1) = 0.815742559275;
      phi_array(2) = 1.15723173726;
      phi_array(3) = 1.39227668949;
      phi_array(4) = 1.58148598466;
      phi_array(5) = 1.74554191992;
      phi_array(6) = 1.89424894726;
      phi_array(7) = 2.03322690864;
      phi_array(8) = 2.16616784027;
      phi_array(9) = 2.29581001175;
      phi_array(10) = 2.42444029861;
      phi_array(11) = 2.55421586793;
      phi_array(12) = 2.68743963225;
      phi_array(13) = 2.82689097563;
      phi_array(14) = 2.97634882193;
      phi_array(15) = 3.14159265359;
      phi_array(16) = 3.33283428007;
      phi_array(17) = 3.57189386259;
      phi_array(18) = 3.92463875127;
      phi_array(19) = 4.80360147994;
      phi_array(20) = 2*PI;
    }

  else
    {
      cerr << "ERROR: Bad index value in getPhiArray." << endl;
      exit(1);
    }

  return phi_array;
}

Array1D<double> AngleDistributions::getThetaArray(int index)
{
  Array1D<double> th_array(21);

  if (index == 1)
    {
      th_array(0) = 0.0;
      th_array(1) = 0.26112746021;
      th_array(2) = 0.373700868288;
      th_array(3) = 0.463397539852;
      th_array(4) = 0.54207933412;
      th_array(5) = 0.614395430249;
      th_array(6) = 0.682799390386;
      th_array(7) = 0.748840434695;
      th_array(8) = 0.813628767581;
      th_array(9) = 0.878063821319;
      th_array(10) = 0.942953024028;
      th_array(11) = 1.00910671532;
      th_array(12) = 1.07741148075;
      th_array(13) = 1.14893652963;
      th_array(14) = 1.22508554452;
      th_array(15) = 1.30785910828;
      th_array(16) = 1.40036860432;
      th_array(17) = 1.50809742635;
      th_array(18) = 1.64255076655;
      th_array(19) = 1.8371410273;
      th_array(20) = PI;
    }

  else if (index == 2)
    {
      th_array(0) = 0.0;
      th_array(1) = 1.30445162629;
      th_array(2) = 1.49904188704;
      th_array(3) = 1.63349522724;
      th_array(4) = 1.74122404927;
      th_array(5) = 1.83373354531;
      th_array(6) = 1.91650710907;
      th_array(7) = 1.99265612396;
      th_array(8) = 2.06418117284;
      th_array(9) = 2.13248593827;
      th_array(10) = 2.19863962956;
      th_array(11) = 2.26352883227;
      th_array(12) = 2.32796388601;
      th_array(13) = 2.3927522189;
      th_array(14) = 2.4587932632;
      th_array(15) = 2.52719722334;
      th_array(16) = 2.59951331947;
      th_array(17) = 2.67819511374;
      th_array(18) = 2.7678917853;
      th_array(19) = 2.88046519338;
      th_array(20) = PI;
    }

  else if (index == 3)
    {
      th_array(0) = 0.0;
      th_array(1) = 0.506913286035;
      th_array(2) = 0.70769274906;
      th_array(3) = 0.859005293588;
      th_array(4) = 0.98600496446;
      th_array(5) = 1.09852964564;
      th_array(6) = 1.20169739354;
      th_array(7) = 1.2986211989;
      th_array(8) = 1.3914258519;
      th_array(9) = 1.48171097387;
      th_array(10) = 1.57079632679;
      th_array(11) = 1.65988167972;
      th_array(12) = 1.75016680169;
      th_array(13) = 1.84297145469;
      th_array(14) = 1.93989526005;
      th_array(15) = 2.04306300795;
      th_array(16) = 2.15558768913;
      th_array(17) = 2.28258736;
      th_array(18) = 2.43389990453;
      th_array(19) = 2.63467936756;
      th_array(20) = PI;
    }

  else if (index == 5)
    {
      th_array(0) = 0.0;
      th_array(1) = 0.241768287031;
      th_array(2) = 0.345875378793;
      th_array(3) = 0.428709148315;
      th_array(4) = 0.501285624277;
      th_array(5) = 0.567896992379;
      th_array(6) = 0.630814517004;
      th_array(7) = 0.691462398921;
      th_array(8) = 0.750855548584;
      th_array(9) = 0.809814534384;
      th_array(10) = 0.869067640752;
      th_array(11) = 0.929334557653;
      th_array(12) = 0.991405886537;
      th_array(13) = 1.0562181404;
      th_array(14) = 1.12499738806;
      th_array(15) = 1.19948064362;
      th_array(16) = 1.28236993143;
      th_array(17) = 1.37842926139;
      th_array(18) = 1.49767301427;
      th_array(19) = 1.66997043121;
      th_array(20) = PI;
    }

  else if (index == 6)
    {
      th_array(0) = 0.0;
      th_array(1) = 0.902754351571;
      th_array(2) = 1.11496045417;
      th_array(3) = 1.26573488799;
      th_array(4) = 1.38892685317;
      th_array(5) = 1.49634039476;
      th_array(6) = 1.59369223076;
      th_array(7) = 1.6842739961;
      th_array(8) = 1.77021220092;
      th_array(9) = 1.85303787533;
      th_array(10) = 1.93393724916;
      th_array(11) = 2.01392608495;
      th_array(12) = 2.09395212875;
      th_array(13) = 2.17499281701;
      th_array(14) = 2.25816573717;
      th_array(15) = 2.34488114462;
      th_array(16) = 2.43713682858;
      th_array(17) = 2.53813250998;
      th_array(18) = 2.65396231898;
      th_array(19) = 2.8001870172;
      th_array(20) = PI;
    }      

  else if (index == 8)
    {
      th_array(0) = 0.0;
      th_array(1) = 0.547750423693;
      th_array(2) = 0.749461867664;
      th_array(3) = 0.897213238891;
      th_array(4) = 1.01927957278;
      th_array(5) = 1.12640029714;
      th_array(6) = 1.22400817797;
      th_array(7) = 1.31533549507;
      th_array(8) = 1.40255342244;
      th_array(9) = 1.48726842322;
      th_array(10) = 1.57079632679;
      th_array(11) = 1.65432423037;
      th_array(12) = 1.73903923115;
      th_array(13) = 1.82625715852;
      th_array(14) = 1.91758447562;
      th_array(15) = 2.01519235645;
      th_array(16) = 2.12231308081;
      th_array(17) = 2.2443794147;
      th_array(18) = 2.39213078593;
      th_array(19) = 2.5938422299;
      th_array(20) = PI;
    }

  else if (index == 11)
    {
      th_array(0) = 0.0;
      th_array(1) = 0.319600760695;
      th_array(2) = 0.457024183973;
      th_array(3) = 0.56626548907;
      th_array(4) = 0.661849495428;
      th_array(5) = 0.749466155368;
      th_array(6) = 0.832103903513;
      th_array(7) = 0.911637753426;
      th_array(8) = 0.989397076186;
      th_array(9) = 1.06644507178;
      th_array(10) = 1.14371784491;
      th_array(11) = 1.222135812;
      th_array(12) = 1.30268676671;
      th_array(13) = 1.38654232073;
      th_array(14) = 1.47520838786;
      th_array(15) = 1.57079632679;
      th_array(16) = 1.67656955891;
      th_array(17) = 1.7981577729;
      th_array(18) = 1.94716690535;
      th_array(19) = 2.15650911998;
      th_array(20) = PI;
    }

  else if (index == 13)
    {
      th_array(0) = 0.0;
      th_array(1) = 0.98508353361;
      th_array(2) = 1.19442574824;
      th_array(3) = 1.34343488069;
      th_array(4) = 1.46502309468;
      th_array(5) = 1.57079632679;
      th_array(6) = 1.66638426573;
      th_array(7) = 1.75505033286;
      th_array(8) = 1.83890588688;
      th_array(9) = 1.91945684159;
      th_array(10) = 1.99787480868;
      th_array(11) = 2.07514758181;
      th_array(12) = 2.1521955774;
      th_array(13) = 2.22995490016;
      th_array(14) = 2.30948875007;
      th_array(15) = 2.39212649823;
      th_array(16) = 2.47974315816;
      th_array(17) = 2.57532716452;
      th_array(18) = 2.68456846962;
      th_array(19) = 2.82199189289;
      th_array(20) = PI;
    }
  
  else
    {
      cerr << "ERROR: Bad index in getThetaArray." << endl;
      exit(1);
    }

  return th_array;
}

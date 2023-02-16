#include <math.h>

// std::vector<double> matMultp(std::vector<std::vector<double>> mat1, std::vector<double> mat2) {
//     size_t r3 = mat1.size();
//     size_t c3 = mat1[0].size();
//     size_t m3 = mat2.size();
//     if (m3!=c3) {pcout << "Matrix dimension does not match (matMultp)\n"; exit(EXIT_FAILURE);}

//     std::vector<double> mat3(r3);
//     for (int i = 0; i < r3; ++i) {
//         mat3[i] = 0;
//         for (int j = 0; j < m3; ++j) {
           

//                 mat3[i] += mat1[i][j] * mat2[j];
           
//         }
//     }
//     return mat3;
// }

//3d version
std::vector<double> matMultp(std::vector<std::vector<double>> mat1, std::vector<double> mat2) {
    size_t r3 = mat1.size();
    size_t c3 = mat1[0].size();
    size_t m3 = mat2.size();
    if (m3!=c3) {pcout << "Matrix dimension does not match (matMultp)\n"; exit(EXIT_FAILURE);}

    std::vector<double> mat3(r3);
for (int i = 0; i < r3; ++i) {
    mat3[i] = 0;
    for (int j = 0; j < c3; ++j) {
        for (int k = 0; k < m3; ++k) {
            mat3[i] += mat1[i][j] * mat2[k];
        }
    }
 } 
 return mat3;
}


// double vecMultp(std::vector<double> mat1, std::vector<double> mat2) {
//     size_t r3 = mat1.size();
//     if (r3!=mat2.size()) {pcout << "Matrix dimension does not match (vecMultp)\n"; exit(EXIT_FAILURE);}

//     double d3 = 0;
//     for (int i = 0; i < r3; ++i) {
//         d3 += mat1[i] * mat2[i];
//     }

//     return d3;
// }

//3d version
double vecMultp(std::vector<std::vector<double>> mat1, std::vector<std::vector<double>> mat2) {
    size_t r3 = mat1.size();
    if (r3!=mat2.size()) {pcout << "Matrix dimension does not match (vecMultp)\n"; exit(EXIT_FAILURE);}

    double d3 = 0;
    for (int i = 0; i < r3; ++i) {
        for (int j = 0; j < mat2[0].size(); ++j) {
            for (int k = 0; k < mat1[0].size(); ++k) {
                d3 += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return d3;
}

// std::vector<double> matSum(std::vector<double> mat1, std::vector<double> mat2) {
//     size_t r3 = mat1.size();
//     if (r3!=mat2.size()) {pcout << "Matrix dimension does not match (matSum)\n"; exit(EXIT_FAILURE);}

//     std::vector<double> mat3(r3);
//     for (int i = 0; i < r3; ++i) {
//         mat3[i] = mat1[i] + mat2[i];
//     }
//     return mat3;
// }

//3d version
std::vector<std::vector<std::vector<double>>> matSum(std::vector<std::vector<std::vector<double>>> mat1,std::vector<std::vector<std::vector<double>>> mat2) {
size_t r1 = mat1.size();
size_t r2 = mat1[0].size();
size_t r3 = mat1[0][0].size();
if (r1!=mat2.size() || r2!=mat2[0].size() || r3!=mat2[0][0].size()) {pcout << "Matrix dimension does not match (matSum)\n";exit(EXIT_FAILURE);}

    std::vector<std::vector<std::vector<double>>> mat3(r1,std::vector<std::vector<double>>(r2,std::vector<double>(r3)));
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < r2; ++j) {
            for (int k = 0; k < r3; ++k) {
                mat3[i][j][k] = mat1[i][j][k] + mat2[i][j][k];
            }
        }
    }
    return mat3;
}


// std::vector<double> tansig(std::vector<double> n) {
//     size_t row = n.size();
//     std::vector<double> mat3(row);

//     for (int i = 0; i < row; ++i) {
//         mat3[i] = 2 / (1 + exp(-2*n[i])) - 1;
//     }
//     return mat3;
// }


//3d version
std::vector<std::vector<std::vector<double>>> tansig(std::vector<std::vector<std::vector<double>>> n) {
    size_t depth = n.size();
    size_t rows = n[0].size();
    size_t cols = n[0][0].size();
    std::vector<std::vector<std::vector<double>>> mat3(depth, std::vector<std::vector<double>>(rows, std::vector<double>(cols)));

    for (size_t k = 0; k < depth; ++k) {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                mat3[k][i][j] = 2 / (1 + exp(-2*n[k][i][j])) - 1;
            }
        }
    }
    return mat3;
}

// void defineNeuralNetwork(std::vector<std::vector<double>> Fin, std::vector<std::vector<double>> &Fout, std::vector<double> &bioR )
// {

//     std::vector<double> b1 = {0.67792810869733710621,-0.9886167579067058897,-0.37003642622919635796,1.187986042338397219,1.4928311621218286476,0.74486892316637776101,1.2366922395968444892,0.58038069834857120011,0.69148356138717936847,0.76749740434725033378};
//     std::vector<double> b2 = {-1.1913964702126025319,0.5315756456217934911,0.51812888918392441262,0.78473922822423625156,-0.04780782970920898628,-0.074872203343969892519,0.14691264129258957416,0.28025089037945333237,-0.25030251589138291513,-0.35615209218765614407};
//     std::vector<double> b3 = {0.83095050244120038929,-0.3825593621459533189,-0.30574795305885948959,-0.19276192855578366814,-0.44085952624946872502,0.033006063256286348462,0.041652187406828918015,0.7305377886218068495,0.020819257857229409719,0.4305093214655514311};
//     std::vector<double> b4 = {-0.026055378788094950976,-0.2945281106715050834,0.085090100040697905226,0.29779965595809959611,0.19527693799985568202,0.12917058855744362189,0.48370574551530176599,0.098492259186305872176,0.47705363420121221774,-0.19315483138355993287};
//     double b5 = 0.85641685816650836571;

//     std::vector<std::vector<double>> IW1_1 ={ {-2.4353481594548980205,1.2453128078005795132},
//                                               {3.0447436282826410014,-0.2160224443494934421},
//                                               {3.3652253034249217656,-2.8426260160350889095},
//                                               {-3.1293074164349627964,-0.15144313699638009552},
//                                               {-4.2466279892081377767,0.048067073437338070363},
//                                               {-2.0644520776813846119,0.20579224152417988081},
//                                               {-2.7118141917732496715,-1.1832463926401328713},
//                                               {0.79348921033598163177,0.61294054741632619798},
//                                               {-0.97016886502106802759,-0.44174306611501107378},
//                                               {-2.6909608875897919056,0.41133381228525300877}
//                                             };
//     std::vector<std::vector<double>> LW2_1 = { {-0.13021009621502233067,0.46442254820360506784,-0.16254619072135081947,0.70479133981267150233,-0.80601268077157350866,0.35014543180417939672,0.56774008699689559876,-0.3906759623889550781,-0.898143664504920336,-0.24754261405271288377},
//                                                {2.518750594262290754,-3.3890342415710121848,-0.71120619706686016848,3.8668651434643219744,4.4018552127422019282,2.2613443499182683816,3.9169814443830763828,-0.60161775818946705563,2.7068183968055965494,2.9795796512994212613},
//                                                {0.42557023743950256334,-0.42910835299689670252,0.006623678487498242326,0.32585250139515931078,0.17252820512663838426,-0.33120755333000762022,0.15016415539439734173,0.13336646976018576294,0.23895690885079068355,-0.22793121065459601149},
//                                                {0.20078429524850449628,-0.38157473464206032032,-0.099047076698002681217,0.070102567431937559683,0.35507260532355533478,-0.54461934713816939624,-0.1968711534564925314,0.60250507440798961589,-0.44426905569738106561,0.24315149459187768155},
//                                                {0.31629747518418027674,-0.3118319044532789075,-0.9123417313631246861,-0.2573673873336732032,0.31267570421009849291,-0.13358228571689653719,-0.55898333022648172275,0.10666045246443203731,-0.66973344485640051715,0.24014261400320235929},
//                                                {0.6924282332619320357,-0.078382125152097775755,-3.4840802350832711376,-0.18172819551184554721,0.1515799541514252502,-0.015265525306761567811,-0.32289229015509640641,0.16159579535671905748,-0.66863753416093285598,0.33647512057550205133},
//                                                {-0.25886547528813286245,0.13633458336918899412,0.2947218735301452841,-0.0023282054626124577476,0.054456311793002872002,0.095008845918454917778,0.28959484682779945697,-0.029137712697683795793,0.49594911858347640043,0.039736674169332318607},
//                                                {0.51073247336451210732,-0.14384984072356604701,-0.15639684424700089904,-0.11088610710943468118,-0.022872385465158150131,0.083135507534787433936,-0.21613475120211508851,0.16091041916938766954,-0.18378098674317042138,0.20554819668917792552},
//                                                {0.25237798811669914789,-0.36105550555846621652,0.011381940728288644782,0.56407634030788489365,0.85647370945278011867,-0.14480711335240797899,0.3046788419351742494,0.0013516875590215597994,0.2757582991239274639,0.38411940065103516995},
//                                                {0.33884964013140711492,-0.23090126343401515263,-0.96122609731222552476,-0.33304000614680201453,-0.14417899498289740712,0.33737575440946965255,-0.4533261036858540205,0.053629664722837928903,-0.19685127443443936612,0.093581538729661503662}
//                                              };
//     std::vector<std::vector<double>> LW3_2 = { {-0.82091872069622850994,0.22438251557365621047,0.5265595653322431291,-0.69973974122533599829,0.37042549249387185517,-0.19308018099405999113,-0.10580395159475251832,-0.35336236161202799755,0.1916518290153929327,0.35086123362869919839},
//                                                {0.31091264379149929908,-0.41007743706175431297,-0.48046580580431064167,0.58664915847258081172,0.22131366753307862849,0.15710960104902832457,0.12075542322889883107,-0.028059878601090128963,0.18240355990986434342,0.12440912996266705048},
//                                                {0.96026291928512541585,-0.15062914356935228066,-0.30695289709264289568,0.26381162551885606327,0.38617277157996904302,-0.38170038724703242439,-0.0058172276994241274573,0.012571174361457283439,-0.38188721451308255128,0.28983069504579134223},
//                                                {-0.43721611603198734519,-0.38226614450786755572,-0.064579274030037991938,0.68636814976872972949,0.29061976100115061161,-0.19957639893847745061,-0.11680693987323419181,0.048298821524873372657,-0.27169183799227214493,0.11861785588166218197},
//                                                {-0.35502458178955448309,5.5157577586010413384,0.61302041852483513118,0.42279811612204359905,-0.10291230816136444359,-0.22076224661081123024,0.31129365685323112656,-0.061542142486668921508,0.63168159727335326803,-0.71462674446375573645},
//                                                {-0.72917896299636197899,6.1331129692495380823,0.94029284802146151367,-0.62536401033617916578,-0.22084173697346587417,-0.19512019308629499625,0.11926434833272134273,0.031419645512173362267,0.57682446333386916404,-0.98336946536787828155},
//                                                {0.098646442334866371593,-0.37893382044506546125,-0.34244456237860165793,0.032650961554379520635,0.41080882671301532927,-0.37502061540207848322,-0.42317490376427258081,0.072047219209698989961,0.18093273528051156962,0.043773091105185450711},
//                                                {-0.4841899677432765503,2.1348341673979480682,0.43061476346751970112,-0.90133856813014867626,-0.64908551956644688907,-0.63958696082321175869,0.7022518391949992278,-0.070022874034757501271,0.028640308837575328277,-0.45938347454213213084},
//                                                {-0.64169883917626657777,1.4295258938535240212,0.58995300762194913258,-0.34490238270280276778,-0.057939343936153311909,-0.41764281575526895907,0.3855552452338121272,0.028627464540422968564,0.030846031053558486956,-0.8270502539059545466},
//                                                {-0.073115059806730969827,0.52703172973593737094,-0.43719128364029113953,0.33410566336205554938,0.050351153505058941773,-0.00099089668102476154837,0.24476931026614226483,0.4556475612623314686,-0.24209632021610011376,-0.25351670878630366834}
//                                              };
//     std::vector<std::vector<double>> LW4_3 = { {0.059376900871753603151,0.065211221136606573046,0.1271323327647514434,0.052177527427534793614,-0.091507903248997726764,-0.083387337494353508394,-0.12177727144466253539,0.09500991491455457183,-0.19593155730913847101,-0.088671200217663559418},
//                                                {-0.089769771890704092021,0.76643193439072965223,0.60761302769193525908,0.12116809552328869359,-0.75750278302614015846,-0.6418405634060270204,0.098320409066735614534,-0.33103336379812819956,0.026766804490196069444,0.11232913636935391855},
//                                                {0.037855696090291823808,0.014062954034362624631,-0.17573398872942558313,0.13200451401746191027,-0.078003996544937156954,-0.036978318665778164842,-0.060842383502702331033,0.027339729725764701923,0.039318031935295143231,-0.15346890150974940026},
//                                                {0.57367122054253938401,-0.37859251209054933796,0.14403494026234808789,-0.61829591949860285283,1.4075355971549599055,1.4353399340805503837,-0.51715996925715501664,1.1563131883992683324,-0.12644624692200742699,-0.44830438379016496198},
//                                                {-0.092141549145137904842,0.92593116399977504205,0.16547819777257360974,0.45047699556308307134,-0.25849350700304479789,-2.6775146159887666109,0.013845636158351385531,-0.45942665476506833189,-0.50885720864770955796,0.17826599281715602152},
//                                                {0.98560263600769482117,-0.6068022719561978473,-0.14405888941267674941,-0.37402020425467558118,4.8060884258287162041,4.6669221352427765481,-0.31396656433810349318,2.6511657592973332243,1.7851072292153427057,0.26026441566828961705},
//                                                {0.037965905954689092849,0.078946346238600167977,-0.051274048622614552817,-0.16525330221492048888,0.2445666393614064904,0.91982039311092222977,0.028837994974132648285,0.24866123979819970691,0.22687031011158714788,0.48090175189056399985},
//                                                {0.25266185203144747584,-0.22712153817671429379,-0.15841465368015558712,-0.066760356633889891831,0.089591994516210693433,-0.10791148810566909833,-0.055028367611793319036,-0.15914571815020828183,0.081741948688057744499,-0.10092094147119363978},
//                                                {1.3040893149970238518,-0.55801226086847832697,-0.98044396835311597993,0.23425980957510209035,0.18021286832592975369,0.2874211884096304348,0.40338674326675910686,0.33142284410054684285,-0.089475538864739537215,-0.32480745373934855058},
//                                                {0.0096018183614487231936,0.062714743404617773193,0.027392093261605989646,0.049642895566454105227,0.037226016709552951778,-0.15834473600114601366,-0.18006539194424373007,0.041681844469807542708,-0.065903487355015166749,-0.096443697730371300003}
//                                              };
//     std::vector<double> LW5_4 = {0.082646451299537612711,0.37850905959073699591,0.059437553346362914652,-1.0173724340607672723,-1.1624070560978583266,-0.84840045855708201561,0.50599701115930084683,0.044083723718335611486,-1.3130720140869056589,-0.10527344255602824608};

//     double x_offset0 = 0.000890159834356918;
//     double x_offset1 = 3.51303254819135e-05;
//     double x_gain0 = 0.200059307889652;
//     double x_gain1 = 4.00231590228802;
//     double x_ymin = -1;

//     double y_ymin = -1;
//     double y_gain = 37.9181676275123;
//     double y_offset = 0;

//     std::vector<double> Xp0(2);
//     Xp0[0] = (Fin[0][0]-x_offset0)*x_gain0+x_ymin;
//     Xp0[1] = (Fin[0][1]-x_offset1)*x_gain1+x_ymin;

//     // Layers
//     std::vector<double> a1 = tansig( matSum(b1, matMultp(IW1_1,Xp0)) );
//     std::vector<double> a2 = tansig( matSum(b2, matMultp(LW2_1,a1)) );
//     std::vector<double> a3 = tansig( matSum(b3, matMultp(LW3_2,a2)) );
//     std::vector<double> a4 = tansig( matSum(b4, matMultp(LW4_3,a3)) );
//     double a5 = b5 + vecMultp(LW5_4,a4);

//     // Output
//     double g = (a5-y_ymin)/y_gain+y_offset;
//     if (g < 1e-8) { g = 0; }

//     bioR[0] = g;

// }


//3D Version
void defineNeuralNetwork3D(std::vector<std::vector<std::vector<double>>> Fin, std::vector<std::vector<std::vector<double>>> &Fout, std::vector<double> &bioR )
{
    std::vector<double> b1 = {0.67792810869733710621,-0.9886167579067058897,-0.37003642622919635796,1.187986042338397219,1.4928311621218286476,0.74486892316637776101,1.2366922395968444892,0.58038069834857120011,0.69148356138717936847,0.76749740434725033378};
    std::vector<double> b2 = {-1.1913964702126025319,0.5315756456217934911,0.51812888918392441262,0.78473922822423625156,-0.04780782970920898628,-0.074872203343969892519,0.14691264129258957416,0.28025089037945333237,-0.25030251589138291513,-0.35615209218765614407};
    std::vector<double> b3 = {0.83095050244120038929,-0.3825593621459533189,-0.30574795305885948959,-0.19276192855578366814,-0.44085952624946872502,0.033006063256286348462,0.041652187406828918015,0.7305377886218068495,0.020819257857229409719,0.4305093214655514311};
    std::vector<double> b4 = {-0.026055378788094950976,-0.2945281106715050834,0.085090100040697905226,0.29779965595809959611,0.19527693799985568202,0.12917058855744362189,0.48370574551530176599,0.098492259186305872176,0.47705363420121221774,-0.19315483138355993287};
    double b5 = 0.85641685816650836571;

    std::vector<std::vector<std::vector<double>>> IW1_1 = { 
                                                    {{-2.4353481594548980205,1.2453128078005795132, 0.0}},
                                                    {{3.0447436282826410014,-0.2160224443494934421, 0.0}},
                                                    {{3.3652253034249217656,-2.8426260160350889095, 0.0}},
                                                    {{-3.1293074164349627964,-0.15144313699638009552, 0.0}},
                                                    {{-4.2466279892081377767,0.048067073437338070363, 0.0}},
                                                    {{-2.0644520776813846119,0.20579224152417988081, 0.0}},
                                                    {{-2.7118141917732496715,-1.1832463926401328713, 0.0}},
                                                    {{0.79348921033598163177,0.61294054741632619798, 0.0}},
                                                    {{-0.97016886502106802759,-0.44174306611501107378, 0.0}},
                                                    {{-2.6909608875897919056,0.41133381228525300877, 0.0}}
                                                    };
    std::vector<std::vector<double>> LW2_1 = { {-0.13021009621502233067,0.46442254820360506784,-0.16254619072135081947,0.70479133981267150233,-0.80601268077157350866,0.35014543180417939672,0.56774008699689559876,-0.3906759623889550781,-0.898143664504920336,-0.24754261405271288377},
                                                {2.518750594262290754,-3.3890342415710121848,-0.71120619706686016848,3.8668651434643219744,4.4018552127422019282,2.2613443499182683816,3.9169814443830763828,-0.60161775818946705563,2.7068183968055965494,2.9795796512994212613},
                                                {0.42557023743950256334,-0.42910835299689670252,0.006623678487498242326,0.32585250139515931078,0.17252820512663838426,-0.33120755333000762022,0.15016415539439734173,0.13336646976018576294,0.23895690885079068355,-0.22793121065459601149},
                                                {0.20078429524850449628,-0.38157473464206032032,-0.099047076698002681217,0.070102567431937559683,0.35507260532355533478,-0.54461934713816939624,-0.1968711534564925314,0.60250507440798961589,-0.44426905569738106561,0.24315149459187768155},
                                                {0.31629747518418027674,-0.3118319044532789075,-0.9123417313631246861,-0.2573673873336732032,0.31267570421009849291,-0.13358228571689653719,-0.55898333022648172275,0.10666045246443203731,-0.66973344485640051715,0.24014261400320235929},
                                                {0.6924282332619320357,-0.078382125152097775755,-3.4840802350832711376,-0.18172819551184554721,0.1515799541514252502,-0.015265525306761567811,-0.32289229015509640641,0.16159579535671905748,-0.66863753416093285598,0.33647512057550205133},
                                                {-0.25886547528813286245,0.13633458336918899412,0.2947218735301452841,-0.0023282054626124577476,0.054456311793002872002,0.095008845918454917778,0.28959484682779945697,-0.029137712697683795793,0.49594911858347640043,0.039736674169332318607},
                                                {0.51073247336451210732,-0.14384984072356604701,-0.15639684424700089904,-0.11088610710943468118,-0.022872385465158150131,0.083135507534787433936,-0.21613475120211508851,0.16091041916938766954,-0.18378098674317042138,0.20554819668917792552},
                                                {0.25237798811669914789,-0.36105550555846621652,0.011381940728288644782,0.56407634030788489365,0.85647370945278011867,-0.14480711335240797899,0.3046788419351742494,0.0013516875590215597994,0.2757582991239274639,0.38411940065103516995},
                                                {0.33884964013140711492,-0.23090126343401515263,-0.96122609731222552476,-0.33304000614680201453,-0.14417899498289740712,0.33737575440946965255,-0.4533261036858540205,0.053629664722837928903,-0.19685127443443936612,0.093581538729661503662}
                                                };
    std::vector<std::vector<std::vector<double>>> LW2_2 = { {-0.82091872069622850994,0.22438251557365621047,0.5265595653322431291,-0.69973974122533599829,0.37042549249387185517,-0.19308018099405999113,-0.10580395159475251832,-0.35336236161202799755,0.1916518290153929327,0.35086123362869919839},
                                                {0.31091264379149929908,-0.41007743706175431297,-0.48046580580431064167,0.58664915847258081172,0.22131366753307862849,0.15710960104902832457,0.12075542322889883107,-0.028059878601090128963,0.18240355990986434342,0.12440912996266705048},
                                                {0.96026291928512541585,-0.15062914356935228066,-0.30695289709264289568,0.26381162551885606327,0.38617277157996904302,-0.38170038724703242439,-0.0058172276994241274573,0.012571174361457283439,-0.38188721451308255128,0.28983069504579134223},
                                                {-0.43721611603198734519,-0.38226614450786755572,-0.064579274030037991938,0.68636814976872972949,0.29061976100115061161,-0.19957639893847745061,-0.11680693987323419181,0.048298821524873372657,-0.27169183799227214493,0.11861785588166218197},
                                                {-0.35502458178955448309,5.5157577586010413384,0.61302041852483513118,0.42279811612204359905,-0.10291230816136444359,-0.22076224661081123024,0.31129365685323112656,-0.061542142486668921508,0.63168159727335326803,-0.71462674446375573645},
                                                {-0.72917896299636197899,6.1331129692495380823,0.94029284802146151367,-0.62536401033617916578,-0.22084173697346587417,-0.19512019308629499625,0.11926434833272134273,0.031419645512173362267,0.57682446333386916404,-0.98336946536787828155},
                                                {0.098646442334866371593,-0.37893382044506546125,-0.34244456237860165793,0.032650961554379520635,0.41080882671301532927,-0.37502061540207848322,-0.42317490376427258081,0.072047219209698989961,0.18093273528051156962,0.043773091105185450711},
                                                {-0.4841899677432765503,2.1348341673979480682,0.43061476346751970112,-0.90133856813014867626,-0.64908551956644688907,-0.63958696082321175869,0.7022518391949992278,-0.070022874034757501271,0.028640308837575328277,-0.45938347454213213084},
                                                {-0.64169883917626657777,1.4295258938535240212,0.58995300762194913258,-0.34490238270280276778,-0.057939343936153311909,-0.41764281575526895907,0.3855552452338121272,0.028627464540422968564,0.030846031053558486956,-0.8270502539059545466},
                                                {-0.073115059806730969827,0.52703172973593737094,-0.43719128364029113953,0.33410566336205554938,0.050351153505058941773,-0.00099089668102476154837,0.24476931026614226483,0.4556475612623314686,-0.24209632021610011376,-0.25351670878630366834}
                                                };
    std::vector<std::vector<std::vector<double>>> LW4_3 = { {0.059376900871753603151,0.065211221136606573046,0.1271323327647514434,0.052177527427534793614,-0.091507903248997726764,-0.083387337494353508394,-0.12177727144466253539,0.09500991491455457183,-0.19593155730913847101,-0.088671200217663559418},
                                                {-0.089769771890704092021,0.76643193439072965223,0.60761302769193525908,0.12116809552328869359,-0.75750278302614015846,-0.6418405634060270204,0.098320409066735614534,-0.33103336379812819956,0.026766804490196069444,0.11232913636935391855},
                                                {0.037855696090291823808,0.014062954034362624631,-0.17573398872942558313,0.13200451401746191027,-0.078003996544937156954,-0.036978318665778164842,-0.060842383502702331033,0.027339729725764701923,0.039318031935295143231,-0.15346890150974940026},
                                                {0.57367122054253938401,-0.37859251209054933796,0.14403494026234808789,-0.61829591949860285283,1.4075355971549599055,1.4353399340805503837,-0.51715996925715501664,1.1563131883992683324,-0.12644624692200742699,-0.44830438379016496198},
                                                {-0.092141549145137904842,0.92593116399977504205,0.16547819777257360974,0.45047699556308307134,-0.25849350700304479789,-2.6775146159887666109,0.013845636158351385531,-0.45942665476506833189,-0.50885720864770955796,0.17826599281715602152},
                                                {0.98560263600769482117,-0.6068022719561978473,-0.14405888941267674941,-0.37402020425467558118,4.8060884258287162041,4.6669221352427765481,-0.31396656433810349318,2.6511657592973332243,1.7851072292153427057,0.26026441566828961705},
                                                {0.037965905954689092849,0.078946346238600167977,-0.051274048622614552817,-0.16525330221492048888,0.2445666393614064904,0.91982039311092222977,0.028837994974132648285,0.24866123979819970691,0.22687031011158714788,0.48090175189056399985},
                                                {0.25266185203144747584,-0.22712153817671429379,-0.15841465368015558712,-0.066760356633889891831,0.089591994516210693433,-0.10791148810566909833,-0.055028367611793319036,-0.15914571815020828183,0.081741948688057744499,-0.10092094147119363978},
                                                {1.3040893149970238518,-0.55801226086847832697,-0.98044396835311597993,0.23425980957510209035,0.18021286832592975369,0.2874211884096304348,0.40338674326675910686,0.33142284410054684285,-0.089475538864739537215,-0.32480745373934855058},
                                                {0.0096018183614487231936,0.062714743404617773193,0.027392093261605989646,0.049642895566454105227,0.037226016709552951778,-0.15834473600114601366,-0.18006539194424373007,0.041681844469807542708,-0.065903487355015166749,-0.096443697730371300003}
                                                };
    std::vector<std::vector<double>> LW5_4 = { 0.082646451299537612711,0.37850905959073699591,0.059437553346362914652,-1.0173724340607672723,-1.1624070560978583266,-0.84840045855708201561,0.50599701115930084683,0.044083723718335611486,-1.3130720140869056589,-0.10527344255602824608};

    double x_offset0 = 0.000890159834356918;
    double x_offset1 = 3.51303254819135e-05;
    double x_offset2 = 0.000182730138990035;
    double x_Gain0 = 0.200059307889652;
    double x_Gain1 = 4.00231590228802;
    double x_Gain2 = 0.300019894767784;
    double x_ymin = -1;

    double y_offset0 = 0.000890159834356918;
    double y_offset1 = 3.51303254819135e-05;
    double y_offset2 = 0.000182730138990035;
    double y_Gain0 = 0.200059307889652;
    double y_Gain1 = 4.00231590228802;
    double y_Gain2 = 0.300019894767784;
    double y_ymin = -1;

    double z_offset0 = 0.000890159834356918;
    double z_offset1 = 3.51303254819135e-05;
    double z_offset2 = 0.000182730138990035;
    double z_Gain0 = 0.200059307889652;
    double z_Gain1 = 4.00231590228802;
    double z_Gain2 = 0.300019894767784;
    double z_ymin = -1;

    std::vector<std::vector<double>> Xp(3, std::vector<double>(3));
    Xp[0][0] = (Fin[0][0]-x_offset0)*x_Gain0+x_ymin;
    Xp[0][1] = (Fin[0][1]-x_offset1)*x_Gain1+x_ymin;
    Xp[0][2] = (Fin[0][2]-x_offset2)*x_Gain2+x_ymin;
    Xp[1][0] = (Fin[1][0]-y_offset0)*y_Gain0+y_ymin;
    Xp[1][1] = (Fin[1][1]-y_offset1)*y_Gain1+y_ymin;
    Xp[1][2] = (Fin[1][2]-y_offset2)*y_Gain2+y_ymin;
    Xp[2][0] = (Fin[2][0]-z_offset0)*z_Gain0+z_ymin;
    Xp[2][1] = (Fin[2][1]-z_offset1)*z_Gain1+z_ymin;
    Xp[2][2] = (Fin[2][2]-z_offset2)*z_Gain2+z_ymin;

    std::vector<double> a1 = tansig( matSum(b1, matMultp(IW1_1,Xp0)) );
    std::vector<double> a2 = tansig( matSum(b2, matMultp(LW2_1,a1)) );
    std::vector<double> a3 = tansig( matSum(b3, matMultp(LW3_2,a2)) );
    std::vector<double> a4 = tansig( matSum(b4, matMultp(LW4_3,a3)) );
    double a5 = b5 + vecMultp(LW5_4,a4);
    double a6 = b6 + vecMultp(LW6_5,a5);
    double a7 = b7 + vecMultp(LW7_5,a5);

    // Output
    double g = (a5-y_ymin)/y_Gain + y_offset;
    double h = (a6-z_ymin)/z_Gain + z_offset;
    double i = (a7-z_ymin)/z_Gain + z_offset;
    if (g < 1e-8) { g = 0; }
    if (h < 1e-8) { h = 0; }
    if (i < 1e-8) { i = 0; }

    bioR[0] = g;
    bioR[1] = h;
    bioR[2] = i;

}

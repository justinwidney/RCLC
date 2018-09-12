#REM" THIS VERSION OF RCLC IS A MODif IED VERSION OF TOM CHACKO'S 1996 ORIGINAL."
#REM" THE MODif ICATIONS HAVE BEEN MADE BY CHRIS MCFARLANE (04-97), "
#REM" DAVID PATTISON (05-97 AND 10-01) AND TOM CHACKO (05-02)."
#REM" "
#REM" RCLC CALCULATES P AND T FROM GRT-PL-OPX-QTZ EQUILIBRIA (GRT-OPX-PL-QTZ BAROMETER"
#REM" AND GRT-OPX AL-IN-OPX THERMOMETER), CORRECTED FOR LATE FE-MG EXCHANGE. "
#REM" "
#REM" FE-MG RATIOS OF GRT, OPX, BT AND CRD ARE ADJUSTED ACCORDING TO THEIR MODES"
#REM" TO SATISFY GRT-OPX, GRT-CRD & GRT-BT FE-MG EXCHANGE EQUILIBRIA."
#REM" "
#REM" ALL THERMODYNAMIC DATA (END-MEMBERS AND SOLUTION MODELS) ARE FROM"
#REM" TWQ202b (FILES BA96A.DAT AND BA96A.SLN), EXCEPT FOR BT WHICH USES"
#REM" THE SOLUTION MODEL OF MCMULLIN (1991, CAN. MINERAL.). FOR REFS, SEE BERMAN (1991,
#REM" CAN. MINERAL.), BERMAN AND ARANOVICH (1996, CONTRIB. MINERAL. PETROL.) AND"
#REM" ARANOVICH AND BERMAN (1997, AM. MINERAL.)"
#REM" "
#REM" REQUIRED INPUT DATA ARE CATIONS IN GARNET (12 OXYGENS), OPX (6 OXY), CRD"
#REM" (18 OXY), BT (11 OXY) AND PL (8 OXY), AND MODAL PROPORTIONS OF GRT, OPX,"
#REM" CRD AND BT. if  THERE IS NO CRD OR BT, ZEROES ARE INSERTED AND THE PROGRAM STILL WORKS."
#REM" if  OPX ANALYSES ARE NOT RECALCULATED FOR FE3+, REPORT TOTAL IRON AS FE2+ IN THE INPUT FILE."


import math
import sys
import csv

def readFile(filename):
    global FEGAR, MNGAR, MGGAR, CAGAR, MODEGAR
    global SIOPX, TIOPX, ALOPX, CROPX, Fe3OPX, FE2OPX, MNOPX, MGOPX, CAOPX, MODEOPX
    global FECRD, MNCRD, MGCRD, MODECRD
    global SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT, MODEBT
    global CA, NA, KPL

    file = open(filename, "r")
    lines = file.read().splitlines()

    dataList = []
    for line in lines:
        if( lines.index(line) % 2 == 1):
            mylist = line.replace(' ','').split(',')
            mylist.pop(0)
            dataList.append(mylist)

    FEGAR, MNGAR, MGGAR, CAGAR, MODEGAR = float(dataList[0][0]), float(dataList[0][1]), float(dataList[0][2]), float(dataList[0][3]), float(dataList[0][4])
    SIOPX, TIOPX, ALOPX, CROPX, Fe3OPX, FE2OPX, MNOPX, MGOPX, CAOPX, MODEOPX = float(dataList[1][0]), float(dataList[1][1]), float(dataList[1][2]), float(dataList[1][3]), float(dataList[1][4]), float(dataList[1][5]), float(dataList[1][6]), float(dataList[1][7]), float(dataList[1][8]), float(dataList[1][9])
    FECRD, MNCRD, MGCRD, MODECRD = float(dataList[2][0]), float(dataList[2][1]), float(dataList[2][2]), float(dataList[2][3])
    SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT, MODEBT = float(dataList[3][0]), float(dataList[3][1]), float(dataList[3][2]), float(dataList[3][2]), float(dataList[3][3]), float(dataList[3][4]), float(dataList[3][5]), float(dataList[3][6]), float(dataList[3][7])
    CA, NA, KPL = float(dataList[4][0]), float(dataList[4][1]), float(dataList[4][2])


    for item in dataList:
        print(item)



def ALOPX1():
    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE

    XFEOPX = FE2OPX / 2
    XMGOPX = MGOPX / 2
    XALM1 = (ALOPX - (2 - SIOPX)) / 2
    ALCHOICE = "MODEL 1  XAlOpx (1-site Opx) = (Al - (2-Si) )/2"
    ALCHOICE2 = "XFeOpx = Fe2+/2    XMgOpx = Mg/2 "
    return


def ALOPX2():
    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE
    XFEOPX = FE2OPX / 2
    XMGOPX = MGOPX / 2
    XALM1 = (ALOPX / 2) / 2
    ALCHOICE = "MODEL 2  XAlOpx (1-site Opx) = (Al/2) / 2"
    ALCHOICE2 = "XFeOpx = Fe2+/2    XMgOpx = Mg/2 "
    return


def ALOPX3():
    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE

    XFEOPX = FE2OPX / (FE2OPX + MGOPX + MNOPX + CAOPX + (ALOPX / 2))
    XMGOPX = MGOPX / (FE2OPX + MGOPX + MNOPX + CAOPX + (ALOPX / 2))
    XALM1 = (ALOPX / 2) / (FE2OPX + MGOPX + MNOPX + CAOPX + (ALOPX / 2))
    ALCHOICE = "MODEL 3  XAlOpx (1-site Opx) = (Al/2) / sum"
    ALCHOICE2 = "where sum = Fe2+ + Mg + Mn + Ca + (Al/2). XFeOpx = Fe2+/sum"
    return


def ALOPX4():
    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE
    XFEOPX = FE2OPX / 2
    XMGOPX = MGOPX / 2
    XALM1 = (ALOPX - Fe3OPX - CROPX - (2 * TIOPX)) / 4
    ALCHOICE = "MODEL 4  XAlOpx (1-site Opx) =  ((Al - Fe3+ - Cr - (2*Ti))/2) / 2"
    ALCHOICE2 = "XFeOpx = Fe2+/2    XMgOpx = MgOpx/2 "
    return


def CP():
    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE
    global DATASET, HALM, HPY, HGR, HAN, HBQ, HEN, HFS, HALOPX, HPHL, HANN, HCRD, HFECRD
    global SALM, SPY, SGR, SAN, SBQ, SEN, SFS, SALOPX, SPHL, SANN, SCRD, SFECRD
    global TK, P

    # REM"CALCULATES ENTHALPIES OF PHASES AT T USING STANDARD STATE ENTHALPIES AND HEAT CAPACITY"
    # REM"math.expRESSIONS FROM TWQ202B - BA96A.DAT OF BERMAN"
    HALM = DATASET[0][0] + ( (DATASET[0][2] * (TK - 298.15)) + ((2 * DATASET[0][ 3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[0][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[0][5] * ((TK**-2) - (298.15**-2))))
    HPY = DATASET[1][0] + ((DATASET[1][2] * (TK - 298.15)) + ((2 * DATASET[1][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[1][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[1][5] * ((TK**-2) - (298.15**-2))))
    HGR = DATASET[2][0] + ((DATASET[2][2] * (TK - 298.15)) + ((2 * DATASET[2][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[2][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[2][5] * ((TK**-2) - (298.15**-2))))
    HAN = DATASET[3][0] + ((DATASET[3][2] * (TK - 298.15)) + ((2 * DATASET[3][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[3][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[3][5] * ((TK**-2) - (298.15**-2))))
    HBQ = DATASET[4][ 0] + ((DATASET[4][2] * (TK - 298.15)) + ((2 * DATASET[4][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[4][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[4][5] * ((TK**-2) - (298.15**-2))))
    HEN = DATASET[5][ 0] + ((DATASET[5][ 2] * (TK - 298.15)) + ((2 * DATASET[5][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[5][ 4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[5][5] * ((TK**-2) - (298.15**-2))))
    HFS = DATASET[6][ 0] + ((DATASET[6][ 2] * (TK - 298.15)) + ((2 * DATASET[6][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[6][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[6][ 5] * ((TK**-2) - (298.15**-2))))
    HALOPX = DATASET[7][0] + ((DATASET[7][2] * (TK - 298.15)) + ((2 * DATASET[7][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[7][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[7][5] * ((TK**-2) - (298.15**-2))))
    HPHL = DATASET[8][ 0] + ((DATASET[8][ 2] * (TK - 298.15)) + ((2 * DATASET[8][ 3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[8][ 4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[8][ 5] * ((TK**-2) - (298.15**-2))))
    HANN = DATASET[9][0] + ((DATASET[9][2] * (TK - 298.15)) + ((2 * DATASET[9][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[9][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[9][5] * ((TK**-2) - (298.15**-2))))
    HCRD = DATASET[10][0] + ((DATASET[10][2] * (TK - 298.15)) + ((2 * DATASET[10][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[10][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[10][5] * ((TK**-2) - (298.15**-2))))
    HFECRD = DATASET[11][ 0] + ((DATASET[11][ 2] * (TK - 298.15)) + ((2 * DATASET[11][3]) * ((TK ** .5) - (298.15 ** .5))) - (DATASET[11][4] * ((1/TK) - (1/298.15))) - (.5 * DATASET[11][5] * ((TK**-2) - (298.15**-2))))

    # REM"CALCULATES ENTROPIES OF PHASES AT T USING STANDARD STATE ENTROPIES AND HEAT CAPACITY"
    # REM"math.expRESSIONS FROM TWQ202B - BA96A.DAT OF BERMAN"
    #print("ABC")
    #print("{}, {}, {}, {}, {}".format(DATASET[0][1], DATASET[0][2], DATASET[0][3], DATASET[0][4], TK))

    SALM = DATASET[0][1] + ((DATASET[0][ 2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[0][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[0][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[0][5]) * ((TK ** -3) - (298.15 ** -3))))
    SPY = DATASET[1][1] + ((DATASET[1][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[1][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[1][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[1][5]) * ((TK ** -3) - (298.15 ** -3))))
    SGR = DATASET[2][1] + ((DATASET[2][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[2][3]) * ((TK ** -.5) - (298.15 ** -.5))) -  ((.5 * DATASET[2][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[2][5]) * ((TK ** -3) - (298.15 ** -3))))
    SAN = DATASET[3][1] + ((DATASET[3][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[3][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[3][4]) * ((TK**-2) - (298.15**-2))) -  (((1 / 3) * DATASET[3][5]) * ((TK ** -3) - (298.15 ** -3))))
    SBQ = DATASET[4][1] + ((DATASET[4][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[4][ 3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[4][ 4]) * ((TK**-2) - (298.15**-2))) -    (((1 / 3) * DATASET[4][5]) * ((TK ** -3) - (298.15 ** -3))))
    #print("******")
    #print(SBQ)
    #print(TK)
    #print("*****")
    SEN = DATASET[5][1] + ((DATASET[5][2] * ((math.log(TK)) - (math.log(298.15)))) -     ((2 * DATASET[5][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[5][ 4]) * ((TK**-2) - (298.15**-2))) -    (((1 / 3) * DATASET[5][ 5]) * ((TK ** -3) - (298.15 ** -3))))
    SFS = DATASET[6][1] + ((DATASET[6][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[6][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[6][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[6][5]) * ((TK ** -3) - (298.15 ** -3))))
    SALOPX = DATASET[7][1] + ((DATASET[7][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[7][3]) * ((TK ** -.5) - (298.15 ** -.5))) -  ((.5 * DATASET[7][4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[7][5]) * ((TK ** -3) - (298.15 ** -3))))
    SPHL = DATASET[8][1] + ((DATASET[8][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[8][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[8][ 4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[8][5]) * ((TK ** -3) - (298.15 ** -3))))
    SANN = DATASET[9][1] + ((DATASET[9][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[9][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[9][4]) * ((TK**-2) - (298.15**-2))) -  (((1 / 3) * DATASET[9][5]) * ((TK ** -3) - (298.15 ** -3))))
    SCRD = DATASET[10][ 1] + ((DATASET[10][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[10][ 3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[10][ 4]) * ((TK**-2) - (298.15**-2))) - (((1 / 3) * DATASET[10][5]) * ((TK ** -3) - (298.15 ** -3))))
    SFECRD = DATASET[11][1] + ((DATASET[11][2] * ((math.log(TK)) - (math.log(298.15)))) - ((2 * DATASET[11][3]) * ((TK ** -.5) - (298.15 ** -.5))) - ((.5 * DATASET[11][ 4]) * ((TK**-2) - (298.15**-2))) -  (((1 / 3) * DATASET[11][5]) * ((TK ** -3) - (298.15 ** -3))))

    return


def VOLUMEPT():
    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE
    global PBARS

    global VALM, VPY, VGR, VAN, VBQ, VEN, VFS, VALOPX, VPHL, VANN, VCRD, VFECRD
    global TK, P

    # "CALCULATES VOLUMES OF PHASES AT P AND T USING STANDARD STATE VOLUMES AND math.expANSION"
    # "AND COMPRESSIBILITY math.expRESSIONS FROM TWQ202B - BA96A.DAT OF BERMAN"
    VALM = 11.524 * (1 + (.0000185989054 * (TK - 298)) + (7.4711E-09 * ((TK - 298) ** 2)) + (-.000000570324 * PBARS) + (
                4.344E-13 * (PBARS ** 2)))

    VPY = 11.311 * (1 + (0.0000225186544 * (TK - 298)) + (3.7044E-09 * ((TK - 298) ** 2)) + (-.000000576209 * PBARS) + (4.42E-13 * (PBARS ** 2)))
    VGR = 12.538 * (1 + (0.0000189942017 * (TK - 298)) + (7.9756E-09 * ((TK - 298) ** 2)) + (-.0000006539136 * PBARS) + (1.635E-12 * (PBARS ** 2)))
    VAN = 10.075 * (1 + (.0000109181141 * (TK - 298)) + (4.1985E-09 * ((TK - 298) ** 2)) + (-.0000012724268 * PBARS) + (3.1762E-12 * (PBARS ** 2)))
    VBQ = 2.37 * (1 + (0 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000012382672 * PBARS) + (7.0871E-12 * (PBARS ** 2)))
    VEN = 3.133 * (1 + (.0000246558172 * (TK - 298)) + (7.467E-09 * ((TK - 298) ** 2)) + (-.0000007493458 * PBARS) + (4.467E-13 * (PBARS ** 2)))
    VFS = 3.295 * (1 + (.0000314064017 * (TK - 298)) + (8.04E-09 * ((TK - 298) ** 2)) + (-.0000009111044 * PBARS) + (3.034E-13 * (PBARS ** 2)))
    VALOPX = 3.093 * (1 + (.0000246558172 * (TK - 298)) + (7.467E-09 * ((TK - 298) ** 2)) + (-.0000007493458 * PBARS) + (4.467E-13 * (PBARS ** 2)))
    VPHL = 14.971 * (1 + (.0000344473262 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000016969784 * PBARS) + (0 * (PBARS ** 2)))
    VANN = 15.487 * (1 + (.0000344473262 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000016969784 * PBARS) + (0 * (PBARS ** 2)))
    VCRD = 23.311 * (1 + (.0000030028742 * (TK - 298)) + (1.8017E-09 * ((TK - 298) ** 2)) + (-.0000011582515 * PBARS) + (0 * (PBARS ** 2)))
    VFECRD = 23.706 * (1 + (.0000042647431 * (TK - 298)) + (0 * ((TK - 298) ** 2)) + (-.0000016998228 * PBARS) + (0 * (PBARS ** 2)))
    return


# REM "SUBROUTINE BERMAN"
def BERMAN():

    global AGR, APY, AAL, GAMMAGAR, PBARS
    global TK, P
    global XCAGAR, XMGGAR, XFEGAR, XMNGAR

    # "THIS SUBROUTINE CALCULATES GARNET ACTIVITIES FOR CA-FE-MG-MN GARNET WITH THE MODEL"
    # "IN TWQ202B - BA96a.SLN OF BERMAN"

    X1 = XCAGAR
    X2 = XMGGAR
    X3 = XFEGAR
    X4 = XMNGAR

    W112 = (85529) - (TK * 18.79) + (PBARS * .21)
    W122 = 50874.9 - (TK * 18.79) + (PBARS * .02)
    W113 = 24025.5 - (TK * 9.43) + (PBARS * .17)
    W133 = 9876.2 - (TK * 9.43) + (PBARS * .09)
    W223 = 1307.4 + (PBARS * .01)
    W233 = 2092.4 + (PBARS * .06)
    W123 = 86852.8 - (TK * 28.22) + (PBARS * .28)
    W124 = 82759.9 - (TK * 28.79) + (PBARS * .1)
    W134 = 7053.9 + (TK * 30.01) + (PBARS * .13)
    W234 = (6361) + (TK * 29.44) + (PBARS * .04)
    W224 = (14558) - (TK * (10)) + (PBARS * .04)
    W244 = (14558) - (TK * (10)) + (PBARS * .04)
    W344 = (158) + (TK * 35.1) + (PBARS * .04)
    W334 = (-(19952)) + (TK * 43.78) + (PBARS * .04)

    ### FIXME: ###
    W114 = 0
    W144 = 0


    TERM1GR = (W112 * ((2 * X1 * X2) - (2 * (X1 ** 2) * X2))) + (W122 * ((X2 ** 2) - (2 * X1 * (X2 ** 2))))
    TERM2GR = (W113 * ((2 * X1 * X3) - (2 * (X1 ** 2) * X3))) + (W133 * ((X3 ** 2) - (2 * X1 * (X3 ** 2))))
    TERM3GR = (W114 * ((2 * X1 * X4) - (2 * (X1 ** 2) * X4))) + (W144 * ((X4 ** 2) - (2 * X1 * (X4 ** 2))))
    TERM4GR = (W223 * ((-2) * (X2 ** 2) * X3)) + (W233 * ((-2) * X2 * (X3 ** 2)))
    TERM5GR = (W224 * ((-2) * (X2 ** 2) * X4)) + (W244 * ((-2) * X2 * (X4 ** 2)))
    TERM6GR = (W334 * ((-2) * (X3 ** 2) * X4)) + (W344 * ((-2) * X3 * (X4 ** 2)))
    TERM7GR = (W123 * ((X2 * X3) - (2 * X1 * X2 * X3))) + (W124 * ((X2 * X4) - (2 * X1 * X2 * X4)))
    TERM8GR = (W134 * ((X3 * X4) - (2 * X1 * X3 * X4))) + (W234 * ((-2) * X2 * X3 * X4))
    SUMGR = TERM1GR + TERM2GR + TERM3GR + TERM4GR + TERM5GR + TERM6GR + TERM7GR + TERM8GR

    ## ** math.exp **
    GAMMAGR = math.exp(SUMGR / (3 * 8.314 * TK))
    AGR = ((X1 * GAMMAGR) ** 3)
    TERM1PY = (W112 * ((X1 ** 2) - (2 * X2 * (X1 ** 2)))) + (W122 * ((2 * X1 * X2) - (2 * (X2 ** 2) * X1)))
    TERM2PY = (W113 * ((-2) * (X1 ** 2) * X3)) + (W133 * ((-2) * X1 * (X3 ** 2)))
    TERM3PY = (W114 * ((-2) * (X1 ** 2) * X4)) + (W144 * ((-2) * X1 * (X4 ** 2)))
    TERM4PY = (W223 * ((2 * X2 * X3) - (2 * (X2 ** 2) * X3))) + (W233 * ((X3 ** 2) - (2 * X2 * (X3 ** 2))))
    TERM5PY = (W224 * ((2 * X2 * X4) - (2 * (X2 ** 2) * X4))) + (W244 * ((X4 ** 2) - (2 * X2 * (X4 ** 2))))
    TERM6PY = (W334 * ((-2) * (X3 ** 2) * X4)) + (W344 * ((-2) * X3 * (X4 ** 2)))
    TERM7PY = (W123 * ((X1 * X3) - (2 * X1 * X2 * X3))) + (W124 * ((X1 * X4) - (2 * X1 * X2 * X4)))
    TERM8PY = (W134 * ((-2) * X1 * X3 * X4)) + (W234 * ((X3 * X4) - (2 * X2 * X3 * X4)))
    SUMPY = TERM1PY + TERM2PY + TERM3PY + TERM4PY + TERM5PY + TERM6PY + TERM7PY + TERM8PY
    GAMMAPY = math.exp(SUMPY / (3 * 8.314 * TK))
    APY = ((X2 * GAMMAPY) ** 3)
    TERM1AL = (W112 * ((-2) * (X1 ** 2) * X2)) + (W122 * ((-2) * X1 * (X2 ** 2)))
    TERM2AL = (W113 * ((X1 ** 2) - (2 * X3 * (X1 ** 2)))) + (W133 * ((2 * X1 * X3) - (2 * (X3 ** 2) * X1)))
    TERM3AL = (W114 * ((-2) * (X1 ** 2) * X4)) + (W144 * ((-2) * X1 * (X4 ** 2)))
    TERM4AL = (W223 * ((X2 ** 2) - (2 * X3 * (X2 ** 2)))) + (W233 * ((2 * X2 * X3) - (2 * (X3 ** 2) * X2)))
    TERM5AL = (W224 * ((-2) * (X2 ** 2) * X4)) + (W244 * ((-2) * X2 * (X4 ** 2)))
    TERM6AL = (W334 * ((2 * X3 * X4) - (2 * (X3 ** 2) * X4))) + (W344 * ((X4 ** 2) - (2 * X3 * (X4 ** 2))))
    TERM7AL = (W123 * ((X1 * X2) - (2 * X1 * X2 * X3))) + (W124 * ((-2) * X1 * X2 * X4))
    TERM8AL = (W134 * ((X1 * X4) - (2 * X1 * X3 * X4))) + (W234 * ((X2 * X4) - (2 * X2 * X3 * X4)))
    SUMAL = TERM1AL + TERM2AL + TERM3AL + TERM4AL + TERM5AL + TERM6AL + TERM7AL + TERM8AL
    GAMMAAL = math.exp(SUMAL / (3 * 8.314 * TK))
    AAL = ((X3 * GAMMAAL) ** 3)
    GAMMAGAR = GAMMAAL / GAMMAPY
    return


# "SUBROUTINE FUHRMAN"
def FUHRMAN():

    global AAN, XAB, XORT, XAN
    global TK, P

    # "THIS SUBROUTINE CALCULATES PLAGIOCLASE ACTIVITIES WITH THE MODEL"
    # "IN TWQ202B - BA96a.SLN OF BERMAN. IT IS THE MODEL OF FUHRMAN AND LINDSLEY (1988)"
    # "WITH WORABAN MODif IED BY BERMAN FOR TWQ202B."
    WABOR = 18.81 - (TK * .0103) + (P * .39)
    WORAB = 27.32 - (TK * .0103) + (P * .39)
    WABAN = 28.226
    WANAB = 8.471
    WANOR = 52.468 - (P * .12)
    WORAN = 47.396
    WORABAN = 100.0455 - (TK * .0103) - (P * .76)
    FIRSTTERM = WORAB * (XAB * XORT * (.5 - XAN - (2 * XAB)))
    SECONDTERM = WABOR * (XAB * XORT * (.5 - XAN - (2 * XORT)))
    THIRDTERM = WORAN * ((2 * XORT * XAN * (1 - XAN)) + (XAB * XORT * (.5 - XAN)))
    FOURTHTERM = WANOR * (((XORT ** 2) * (1 - (2 * XAN))) + (XAB * XORT * (.5 - XAN)))
    FifTHTERM = WABAN * ((2 * XAB * XAN * (1 - XAN)) + (XAB * XORT * (.5 - XAN)))
    SIXTHTERM = WANAB * (((XAB ** 2) * (1 - (2 * XAN))) + (XAB * XORT * (.5 - XAN)))
    SEVENTHTERM = WORABAN * (XORT * XAB * (1 - (2 * XAN)))
    AAN = math.exp((FIRSTTERM + SECONDTERM + THIRDTERM + FOURTHTERM + FifTHTERM + SIXTHTERM + SEVENTHTERM) / (.008314 * TK))
    AAN = (XAN * (((1 + XAN) ** 2) / 4)) * AAN
    return


# "SUBROUTINE MCMULLIN"
def MCMULLIN():

    global GAMMABT
    global TK, XFEBT, XALBT, XTIBT, XMGBT

    # "THIS SUBROUTINE CALCULATES ANNITE AND PHLOGOPITE ACTIVITIES WITH THE MODEL"
    # "OF MCMULLIN (1991). IT IS THUS Dif FERENT FROM THE MODEL IN TWQ202B - BA96a.SLN OF BERMAN."
    WMGFE = 0
    WMGTI = 58.865
    WMGAL = 75
    WFETI = 30.921
    WFEAL = 63.721
    WTIAL = 0
    RTGAMMAMGBT = ((XFEBT ** 2) * WMGFE) + ((XTIBT ** 2) * WMGTI) + ((XALBT ** 2) * WMGAL) + (XFEBT * XTIBT * (WMGFE + WMGTI - WFETI)) + (XFEBT * XALBT * (WMGFE + WMGAL - WFEAL)) + (XTIBT * XALBT * (WMGTI + WMGAL - WTIAL))
    RTGAMMAFEBT = ((XMGBT ** 2) * WMGFE) + ((XTIBT ** 2) * WFETI) + ((XALBT ** 2) * WFEAL) + (XMGBT * XTIBT * (WMGFE + WFETI - WMGTI)) + (XMGBT * XALBT * (WMGFE + WFEAL - WMGAL)) + (XTIBT * XALBT * (WFETI + WFEAL - WTIAL))
    GAMMAMGBT = math.exp(RTGAMMAMGBT / (.008314 * TK))
    GAMMAFEBT = math.exp(RTGAMMAFEBT / (.008314 * TK))
    GAMMABT = GAMMAMGBT / GAMMAFEBT
    return


# "SUBROUTINE CORDIERITE"
def CORDIERITE():

    global GAMMACRD
    global TK, MGRATIOCRD
    # "THIS SUBROUTINE CALCULATES MGCRD AND FECRD ACTIVITIES WITH THE MODEL"
    # REM"IN TWQ202B - BA96a.SLN OF BERMAN"
    WCRD = -1754.7
    RTGAMMAMGCRD = WCRD * ((1 - MGRATIOCRD) ** 2)
    RTGAMMAFECRD = WCRD * (MGRATIOCRD ** 2)
    GAMMAMGCRD = math.exp(RTGAMMAMGCRD / (8.314 * TK))
    GAMMAFECRD = math.exp(RTGAMMAFECRD / (8.314 * TK))
    GAMMACRD = GAMMAMGCRD / GAMMAFECRD
    return


# "SUBROUTINE ARANOVICH"
def ARANOVICH():

    global AEN, AFS, AALOPX, GAMMAOPX
    global PBARS, TK
    global FERATIOOPX
    global XALM1, XMGOPX, XFEOPX
    # "THIS SUBROUTINE CALCULATES OPX ACTIVITIES WITH THE MODEL"
    # "IN TWQ202B - BA96a.SLN OF BERMAN. IT IS BASED ON THE MODEL OF ARANOVICH AND BERMAN (1997)"
    # "NOTE THAT THE ALOPX ACTIVITY MODEL INCLUDES A DARKEN CORRECTION OF THE FORM"
    # "RTLN(GAMMA)ALOPX =RTLN(GAMMA)ALOPX + FE/(FE+MG)*(DH-T*DS)"

    W12 = -4543.8 + (TK * 3.36)
    W23 = -32213.3 - (PBARS * .69)
    W13 = -26944.5 - (PBARS * .58)

    #-23139.5

    RTGAMMAMGOPX = (W12 * (XFEOPX - (XFEOPX * XMGOPX))) - (W23 * XFEOPX * XALM1) + (W13 * (XALM1 - (XMGOPX * XALM1)))
    RTGAMMAFEOPX = (W12 * (XMGOPX - (XFEOPX * XMGOPX))) + (W23 * (XALM1 - (XFEOPX * XALM1))) - (W13 * XMGOPX * XALM1)
    RTGAMMAALOPX = (-1 * (W12 * XFEOPX * XMGOPX)) + (W23 * (XFEOPX - (XFEOPX * XALM1))) + (
                W13 * (XMGOPX - (XMGOPX * XALM1))) + (FERATIOOPX * ((24307) - (TK * 14.404) + (PBARS * .185)))

    GAMMAMGOPX = math.exp(RTGAMMAMGOPX / (8.314 * TK))
    GAMMAFEOPX = math.exp(RTGAMMAFEOPX / (8.314 * TK))
    #
    # print("----------")
    # print(RTGAMMAALOPX)
    # print(XFEOPX)
    # print(XMGOPX)
    # print(XALM1)
    # print(FERATIOOPX)
    #
    # print("TK=" + str(TK) + " PBARS=" + str(PBARS))
    # print("RTGAMMAALOPX")
    # print(RTGAMMAALOPX)
    #
    # print("GAMMAALOPX")
    # print(math.exp(RTGAMMAALOPX / (8.314 * TK) ))
    # print("----------")




    GAMMAALOPX = math.exp( RTGAMMAALOPX / (8.314 * TK) )


    AEN = XMGOPX * GAMMAMGOPX
    AFS = XFEOPX * GAMMAFEOPX
    AALOPX = XALM1 * GAMMAALOPX
    GAMMAOPX = GAMMAMGOPX / GAMMAFEOPX
    return




# MODEL:
def Model():


    global XMNGAR, XFEGAR, XCAGAR, XMGGAR, XFEGARI, MGRATIOGA, MGRATIOGARI, MODXCAGAR
    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE

    global FEGAR, MNGAR, MGGAR, CAGAR, MODEGAR
    global SIOPX, TIOPX, ALOPX, CROPX, Fe3OPX, FE2OPX, MNOPX, MGOPX, CAOPX, MODEOPX
    global FECRD, MNCRD, MGCRD, MODECRD
    global CA, NA, KPL
    global SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT, MODEBT

    global DATASET, HALM, HPY, HGR, HAN, HBQ, HEN, HFS, HALOPX, HPHL, HANN, HCRD, HFECRD
    global SALM, SPY, SGR, SAN, SBQ, SEN, SFS, SALOPX, SPHL, SANN, SCRD, SFECRD
    global VALM, VPY, VGR, VAN, VBQ, VEN, VFS, VALOPX, VPHL, VANN, VCRD, VFECRD
    global AGR, APY, AAL, GAMMAGAR
    global MGRATIOCRD
    global TK, XFEBT, XALBT, XTIBT, XMGBT

    global FERATIOOPX

    global AAN, XAB, XORT, XAN
    global PBARS, TK, P
    global ANSWER1, TITLE





    print ("")
    print ("YOU HAVE A CHOICE FOR CALCULATING XALM IN OPX.")
    print ("THE FOLLOWING FORMULAE ASSUME A 6-OXYGEN OPX FORMULA.")
    print ("")
    print ("1: XALM = Al - (2 - Si)")
    print ("2: XALM = Al/2 ")
    print ( "3: XALM = (Al/2) / (Fe2+ + Mg + Mn + Ca + (Al/2) )" )
    print ("4: XALM = (Al - Fe3+ - Cr - (2*Ti) ) / 2")
    print ("")

    num = int( input("Please enter 1,2,3,4: "))

    if  num == 1:
        ALOPX1()
    if  num == 2:
        ALOPX2()
    if  num == 3:
        ALOPX3()
    if  num == 4:
        ALOPX4()

    XFEOPXI = XFEOPX
    MGRATIOOPX = MGOPX / (MGOPX + FE2OPX)
    MGRATIOOPXI = MGRATIOOPX
    FERATIOOPX = 1 - MGRATIOOPX
    TOTOPX = XFEOPX + XMGOPX + XALM1+(TIOPX/2)+(Fe3OPX/2)+(CROPX/2)+(MNOPX/2)+(CAOPX/2)

    # Garnet mole fraction calculations
    # Step 2 #
    XMGGAR = MGGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XFEGAR = FEGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XCAGAR = CAGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XMNGAR = MNGAR / (MGGAR + CAGAR + FEGAR + MNGAR)
    XFEGARI = XFEGAR
    MGRATIOGAR = MGGAR / (MGGAR + FEGAR)
    MGRATIOGARI = MGRATIOGAR
    MODXCAGAR = (CAGAR + MNGAR) / (CAGAR + MNGAR + FEGAR + MGGAR)

    # "BIOTITE MOLE FRACTION CALCULATIONS"

    if  (MODEBT < 0.01 ):
            XFEBT = 0
            XMGBT = 0
            XTIBT = 0
            XALBT = 0
            MGRATIOBT = 0
    else:
            ALIVBT = 4.0 - SIBT
            ALVIBT = ALBT - ALIVBT
            XFEBT = FEBT/  (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            XMGBT = MGBT/  (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            XALBT = ALVIBT/ (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            XTIBT = TIBT/  (FEBT+MGBT+ALVIBT+TIBT+MNBT)
            MGRATIOBT = MGBT / (MGBT + FEBT)

    XFEBTI =XFEBT
    MGRATIOBTI = MGRATIOBT


    # "CORDIERITE MOLE FRACTION CALCULATIONS"
    if  (MODECRD < 0.01):
            XFECRD = 0
            XMGCRD = 0
            XMNCRD = 0
            MGRATIOCRD = 0
    else:
            XFECRD = FECRD / (FECRD + MGCRD + MNCRD)
            XMGCRD = MGCRD / (FECRD + MGCRD + MNCRD)
            XMNCRD = MNCRD / (FECRD + MGCRD + MNCRD)
            MGRATIOCRD = MGCRD / (MGCRD + FECRD)

    MGRATIOCRDI = MGRATIOCRD
    XFECRDI = XFECRD

    "PLAGIOCLASE MOLE FRACTIONS"
    XAN = CA / (CA + NA + KPL)
    XAB = NA / (NA + CA + KPL)
    XORT = KPL / (NA + CA + KPL)


    f = open('export.txt','a')



    if  (ANSWER1 != "y" or ANSWER1 != "Y" ):
        f.write("\n")
        f.write("RCLC: GRT-OPX-PL-QTZ P-T ESTIMATES CORRECTED FOR LATE FE-MG EXCHANGE \n")
        f.write("THERMODYNAMIC DATA AND math.expRESSIONS FROM TWQ202B\n")
        f.write("Description: ")
        f.write(TITLE) #print File in and Title
        f.write( "\n")




    # ** Need Clarification **

    if  (ANSWER1 == "y" or ANSWER1 == "Y"):
            f.write(" ")
    if  (ANSWER1 == "y" or ANSWER1 == "Y"):
            f.write(" ")

        # Variable name or Text
    f.write(ALCHOICE)
    f.write(ALCHOICE2)
    f.write("\n")


        #f.write(USING "ModeGrt = ##.## ";MODEGAR; USING " ModeOpx = ##.## ";MODEOPX; USING " ModeBt  = ##.## "; MODEBT; USING " ModeCrd = ##.## "; MODECRD)
    #f.write("ModeGrt = " + str(MODEGAR) + " " + "ModeOpx = " + str(MODEOPX) + " " + "ModeBt = " + str(MODEBT) + " " + " ModeCrd = " + str(MODECRD) + "\n")
    f.write("ModeGrt = {:0.3f}      ModeOpx = {:0.3f}       ModeBt = {:0.3f}        ModeCrd = {:0.3f}    \n".format(MODEGAR, MODEOPX, MODEBT, MODECRD))

        #f.write( USING "XFeGrt  = #.### "; XFEGAR; USING " XMgGrt  = #.### "; XMGGAR; USING " XMnGrt  = #.### "; XMNGAR; USING " XCaGrt  = #.### "; XCAGAR)
    #f.write( "XFeGRT = " + str(XFEGAR) + " " + "XMgGrt = " + str(XMGGAR) + " " + "XMnGrt = " + str(XMNGAR) + " " + "XcaGrt = " + str(XCAGAR) + "\n")
    f.write("XFeGRT = {:0.3f}      XMgGrt = {:0.3f}       XMnGrt = {:0.3f}        XcaGrt = {:0.3f}    \n".format(XFEGAR, XMGGAR, XMNGAR, XCAGAR))

        #f.write(USING "XFeOpx  = #.### "; XFEOPX; USING " XMgOpx  = #.### "; XMGOPX; USING " XAlOpx  = #.### "; XALM1 ; USING " TotOpx  = #.### "; TOTOPX)
    #f.write("XFeOpx = " + str(XFEOPX) + " " + "XMgOpx  =  " + str(XMGOPX) + " " + " XAlOpx =  " + str(XALM1) + " " + " TotOpx  =  " + str(TOTOPX) + "\n")
    f.write("XFeOpx = {:0.3f}     XMgOpx = {:0.3f}        XAlOpx = {:0.3f}        TotOpx = {:0.3f}    \n".format(XFEOPX, XMGOPX, XALM1, TOTOPX))

    if  (MODEBT > 0.01):
            #f.write(USING "XFeBt   = #.### "; XFEBT;  USING " XMgBt   = #.### "; XMGBT ; USING " XAlBt   = #.### "; XALBT ; USING " XTiBt   = #.### "; XTIBT)
        #f.write("XFeBt   ="+ str(XFEBT) +  " " +" XMgBt   =  "+ str(XMGBT) + " " + "XAlBt =" + str(XALBT) + " " + "XTiBt  = " + str(XTIBT) + "\n")
        f.write("XFeBt  = {:0.3f}     XMgBt = {:0.3f}        XAlBt = {:0.3f}        XTiBt = {:0.3f}    \n".format(XFEBT, XMGBT, XALBT, XTIBT))


    if  (MODECRD> 0.01):
            #f.write(USING "XFeCrd  = #.### "; XFECRD; USING " XMgCrd  = #.### "; XMGCRD; USING " XMnCrd  = #.### "; XMNCRD)
        #f.write("XFeCrd  =  " + str(XFECRD) +  " XMgCrd  =  "+ str(XMGCRD) +  " XMnCrd  =  "+ str(XMNCRD) + "\n")
        f.write("XFeCrd  = {:0.3f}    XMgCrd = {:0.3f}       XMnCrd = {:0.3f}    \n".format(XFECRD, XMGCRD, XMNCRD))

        #f.write( USING "XAnPl   = #.### "; XAN;    USING " XAbPl   = #.### "; XAB;)
    #f.write( "XAnPl   =  " + str(XAN) +  " XAbPl   =  " + str(XAB) )
    f.write("XAnPl  = {:0.3f}    XAbPl  = {:0.3f}   \n".format(XAN, XAB))




    # STEP 1 #
    #REM " VOLUME FRACTIONS OF FE-MG MINERALS FROM MODE"
    VFGAR = MODEGAR / (MODEGAR + MODEOPX + MODECRD + MODEBT)
    VFOPX = MODEOPX / (MODEGAR + MODEOPX + MODECRD + MODEBT)
    VFCRD = MODECRD / (MODEGAR + MODEOPX + MODECRD + MODEBT)
    VFBT =  MODEBT  / (MODEGAR + MODEOPX + MODECRD + MODEBT)

    #REM "CONVERT VOLUME FRACTION MINERALS TO MOLE FRACTION MINERALS"

    #REM  "DENSITIES"
    DENSFEGAR = 4.33
    DENSMGGAR = 3.54
    DENSCAGAR = 3.56
    DENSMNGAR = 4.19
    DENSFEOPX = 3.96
    DENSMGOPX = 3.21
    DENSMGCRD = 2.53
    DENSFECRD = 2.78
    DENSFEBT =  3.3
    DENSMGBT =  2.7
    DENSGAR = (DENSFEGAR * XFEGAR) + (DENSMGGAR * XMGGAR) + (DENSCAGAR * XCAGAR) + (DENSMNGAR * XMNGAR)
    DENSOPX = (DENSFEOPX * (1 - MGRATIOOPX)) + (DENSMGOPX * MGRATIOOPX)
    DENSCRD = (DENSFECRD * (1 - MGRATIOCRD)) + (DENSMGCRD * MGRATIOCRD)
    DENSBT =  (DENSFEBT  * (1 - MGRATIOBT))  + (DENSMGBT  * MGRATIOBT)

    #REM  "MOLECULAR WEIGHTS"
    MWOPX = (SIOPX * 28.1) + (TIOPX*47.9) + (ALOPX * 26.1) + (CROPX*(52)) + (Fe3OPX*55.8) + (FE2OPX * 55.8) + (MGOPX * 24.3) + (MNOPX * 54.9) + (CAOPX * 40.1) + (6 * 16)
    MWGAR = (3.00  * 28.1) + (2.00  * 26.1) + (FEGAR * 55.8) + (MGGAR * 24.3) + (MNGAR * 54.9) + (CAGAR * 40.1) + (12 * 16)
    MWCRD = (5.00  * 28.1) + (4.00  * 26.1) + (FECRD * 55.8) + (MGCRD * 24.3) + (MNCRD * 54.9) + (18 * 16)
    MWBT =  (SIBT * 28.1)  + (TIBT * 47.9)  + (ALBT * 26.1)  + (FEBT * 55.8)  + (MNBT * 54.9)  + (MGBT * 24.3)  + (NABT * 23) + (KBT * 39.1) + (11 * 16) + 2

    #REM   "MOLES OF MINERALS AND MOLES OF FE-MG COMPONENTS OF MINERALS"
    MOLEGAR = (VFGAR * DENSGAR) / MWGAR
    MOLEOPX = (VFOPX * DENSOPX) / MWOPX



    if (MODECRD < 0.01):
        MOLECRD = 0
    else:
        MOLECRD = (VFCRD * DENSCRD) / MWCRD
    #END if

    if(MODEBT < 0.01):
        MOLEBT = 0
    else:
        MOLEBT = (VFBT * DENSBT) / MWBT
    #ENDIF




    MOLEFEMGGAR = MOLEGAR * (FEGAR + MGGAR)
    MOLEFEMGOPX = MOLEOPX * (FE2OPX + MGOPX)
    MOLEFEMGCRD = MOLECRD * (FECRD + MGCRD)
    MOLEFEMGBT =  MOLEBT  * (FEBT  + MGBT )

    #REM "MOLE FRACTION OF FE-MG COMPONENTS OF MINERALS"
    MFGAR = MOLEFEMGGAR / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)
    MFOPX = MOLEFEMGOPX / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)
    MFCRD = MOLEFEMGCRD / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)
    MFBT =  MOLEFEMGBT  / (MOLEFEMGGAR + MOLEFEMGOPX + MOLEFEMGCRD + MOLEFEMGBT)

    #REM  "CALCULATE XMG ROCK"
    XMGROCK = (MGRATIOGAR * MFGAR) + (MGRATIOOPX * MFOPX) + (MGRATIOCRD * MFCRD) + (MGRATIOBT * MFBT)
    XMGROCKI = XMGROCK

    #REM  "CALCULATE GRT-OPX FE-MG  -  GRT-OPX-PL-QTZ (FE-END MEMBER)INTERSECTION"
    #REM   "ASSUME INITIAL P-T TO BEGIN. ITERATE 10 TIMES."

    TK = 1000
    PBARS = 3000
    P = 3
    CP()

    #print("S Values")
    #print("SAN = " + str(SAN) )
    #print("SFS = " + str(SFS) )
    #print("SBQ = " + str(SBQ))
    #print("SALM= " + str(SALM))
    #print("SGR= " + str(SGR))

    for J in range(10):
        #REM "CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE"

        # FIXME: Causing an error on numbers,

        BERMAN()
        FUHRMAN()
        ARANOVICH()
        VOLUMEPT()

        #  SAN, SFS, SBQ, SALM, SGR

        DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
        DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
        DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)




        #print("DELTA VALUES")
        #print("DELTAHGAPES = " + str(DELTAHGAPES))
        #print("DELTASGAPES = " + str(DELTASGAPES))
        #print("DELTAVGAPES = " + str(DELTAVGAPES))

        KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))
        #print("KGAPES = " + str(KGAPES))

        P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (math.log(KGAPES)))) / DELTAVGAPES
        PBARS = P * 1000

        #REM "CALCULATE GRT-OPX FE-MG EXCHANGE TEMP AT THIS PRESSURE"
        DELTAHFEMGOPX = (((1 * HEN) + ((1 / 3) * HALM)) - ((1 * HFS) + ((1 / 3) * HPY))) / 1000
        DELTASFEMGOPX = (((1 * SEN) + ((1 / 3) * SALM)) - ((1 * SFS) + ((1 / 3) * SPY))) / 1000
        DELTAVFEMGOPX = ((1 * VEN) + ((1 / 3) * VALM)) - ((1 * VFS) + ((1 / 3) * VPY))
        GAMMAFEMGOPX = GAMMAGAR * GAMMAOPX
        KDGAROPX = (XFEGAR * XMGOPX) / (XMGGAR * XFEOPX)
        TGAROPX = (DELTAHFEMGOPX + (P * DELTAVFEMGOPX)) / (DELTASFEMGOPX - (.008314 * math.log(KDGAROPX)) - (.008314 * math.log(GAMMAFEMGOPX)))
        TK = TGAROPX
        TCGAROPX = TGAROPX - 273
        CP()
        #print("S Values")
        #print("SAN = " + str(SAN) )
        #print("SFS = " + str(SFS) )
        #print("SBQ = " + str(SBQ))
        #print("SALM= " + str(SALM))
        #print("SGR= " + str(SGR))

    # DONE LOOP#

    TGAROPXI = TCGAROPX
    PFEMGGTOPI = P

    if  (MODECRD > 0.01):
        #REM  "CALCULATE GRT-CRD FE-MG  -  GRT-OPX-PL-QTZ (FE END MEMBER) INTERSECTION"
        #REM   "ASSUME INITIAL P-T TO BEGIN. ITERATE 10 TIMES."
        TK = 1000
        PBARS = 3000
        P = 3
        CP()

        for J in range(10):

            #REM "CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE"
            BERMAN()
            FUHRMAN()
            ARANOVICH()
            VOLUMEPT()

            DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
            DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
            DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)
            KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))
            P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (math.log(KGAPES)))) / DELTAVGAPES
            PBARS = P * 1000

            #REM "CALCULATE GRT-CRD FE-MG EXCHANGE TEMP AT THIS PRESSURE"
            CORDIERITE()
            DELTAHFEMGCRD = (((.5 * HCRD) + ((1 / 3) * HALM)) - ((.5 * HFECRD) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGCRD = (((.5 * SCRD) + ((1 / 3) * SALM)) - ((.5 * SFECRD) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGCRD = ((.5 * VCRD) + ((1 / 3) * VALM)) - ((.5 * VFECRD) + ((1 / 3) * VPY))
            GAMMAFEMGCRD = GAMMAGAR * GAMMACRD
            KDGARCRD = (XFEGAR * XMGCRD) / (XMGGAR * XFECRD)
            TGARCRD = (DELTAHFEMGCRD + (P * DELTAVFEMGCRD)) / (DELTASFEMGCRD - (.008314 * math.log(KDGARCRD)) - (.008314 * math.log(GAMMAFEMGCRD)))
            TK = TGARCRD
            TCGARCRD = TGARCRD - 273
            CP()

        TCGARCRDI = TCGARCRD
        PGARCRDI = P

    else:
        TCGARCRD = 0
        TCGARCRDI = 0
        PGARCRD = 0
        PGARCRDI = 0

    #END if


    if  (MODEBT > 0.01):

        #REM  "CALCULATE GRT-BT FE-MG  -  GRT-OPX-PL-QTZ (FE END MEMBER) INTERSECTION"
        #REM   "ASSUME INITIAL P-T TO BEGIN. ITERATE 10 TIMES."

        TK = 1000
        PBARS = 3000
        P = 3
        CP()

        for J in range (10):

            #REM "CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE"
            BERMAN()
            FUHRMAN()
            ARANOVICH()
            VOLUMEPT()
            DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
            DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
            DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)
            KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))

            P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (math.log(KGAPES)))) / DELTAVGAPES
            PBARS = P * 1000

            #REM "CALCULATE GRT-BT FE-MG EXCHANGE TEMP AT THIS PRESSURE"
            MCMULLIN()
            DELTAHFEMGBT = ((((1 / 3) * HPHL) + ((1 / 3) * HALM)) - (((1 / 3) * HANN) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGBT = ((((1 / 3) * SPHL) + ((1 / 3) * SALM)) - (((1 / 3) * SANN) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGBT = (((1 / 3) * VPHL) + ((1 / 3) * VALM)) - (((1 / 3) * VANN) + ((1 / 3) * VPY))
            GAMMAGARBT = GAMMABT * GAMMAGAR
            KDGARBT = (XFEGAR * XMGBT) / (XMGGAR * XFEBT)
            TGARBT = (DELTAHFEMGBT + (P * DELTAVFEMGBT)) / (DELTASFEMGBT - (.008314 * math.log(KDGARBT)) - (.008314 * math.log(GAMMAGARBT)))
            TK = TGARBT
            TCGARBT = TGARBT - 273
            CP()


        TCGARBTI = TCGARBT
        PGARBTI = P

    else:
        TCGARBT = 0
        TCGARBTI = 0
        PGARBT = 0
        PGARBTI = 0
    # END if

    #REM "CALCULATE CONVERGED INTERSECTION OF GRT-OPX AL-SOLUBILITY AND GRT-OPX-PL-QTZ USING"
    #REM "FE-END MEMBER math.expRESSIONS."
    #REM "CONVERGENCE APPROACH - 1. CALCULATE INITIAL INTERSECTION OF GRT-OPX AL-SOLUB"
    #REM "AND GRT-OPX-PL-QTZ."
    #REM "2. CHANGE KD GRT-OPX (AND if  APPLICABLE KD GRT-CRD AND KD GRT-BT) SO COINCIDES WITH 1.
    #REM "3. ADJUST FE/MG RATIOS OF FE-MG MINERALS TO SATif Y KD'S.
    #REM "4. REPEAT 10 TIMES (I = 1 TO 10) TO GET CONVERGENCE.
    #REM "ASSUME INITIAL TEMP TO BEGIN"

    TK = 1000
    PBARS = 3000
    P = 3

    for I in range(10):
        CP()

        #REM "CALCULATE INTERSECTION OF FE-AL-OPX AND GRT-OPX-PL-QTZ IN 10 ITERATIONS (J = 1 TO 10)"

        for J in range (10):

            #REM "CALCULATE GRT-OPX-PL-QTZ (FE-END MEMBER) PRESSURE"
            BERMAN()
            FUHRMAN()
            ARANOVICH()
            VOLUMEPT()
            DELTAHGAPES = (((3 * HAN) + (6 * HFS)) - ((3 * HBQ) + (2 * HALM) + HGR)) / 1000
            DELTASGAPES = (((3 * SAN) + (6 * SFS)) - ((3 * SBQ) + (2 * SALM) + SGR)) / 1000
            DELTAVGAPES = ((3 * VAN) + (6 * VFS)) - ((3 * VBQ) + (2 * VALM) + VGR)
            KGAPES = ((AAN ** 3) * (AFS ** 6)) / (AGR * (AAL ** 2))
            P = ((TK * DELTASGAPES) - DELTAHGAPES - (.008314 * TK * (math.log(KGAPES)))) / DELTAVGAPES
            PBARS = P * 1000

            #REM "CALCULATE FE-ALOPX TEMPERATURE AT THIS PRESSURE"
            DELTAHFEAL = ((HALOPX + (3 * HFS)) - HALM) / 1000
            DELTASFEAL = ((SALOPX + (3 * SFS)) - SALM) / 1000
            DELTAVFEAL = (VALOPX + (3 * VFS)) - VALM
            KFEAL = ((AFS ** 3) * AALOPX) / AAL
            TFEAL = (DELTAHFEAL + (P * DELTAVFEAL)) / (DELTASFEAL - (.008314 * math.log(KFEAL)))
            TK = TFEAL
            TC = TK - 273
            CP()
            #NEXT J


        if  (I==0):
            TFEALI = TC
            PFEALI = P

        #END if

        BERMAN()
        ARANOVICH()
        if  (MODECRD > 0.01):
            CORDIERITE()

        if  (MODEBT > 0.01):
            MCMULLIN()

        VOLUMEPT()

        #REM "CALCULATES A CORRECTED KD(GRT-OPX(FE-MG)) "
        DELTAHFEMGOPX = (((1 * HEN) + ((1 / 3) * HALM)) - ((1 * HFS) + ((1 / 3) * HPY))) / 1000
        DELTASFEMGOPX = (((1 * SEN) + ((1 / 3) * SALM)) - ((1 * SFS) + ((1 / 3) * SPY))) / 1000
        DELTAVFEMGOPX = ((1 * VEN) + ((1 / 3) * VALM)) - ((1 * VFS) + ((1 / 3) * VPY))
        GAMMAFEMGOPX = GAMMAGAR * GAMMAOPX
        KDGAROPX = ((TK * DELTASFEMGOPX) - DELTAHFEMGOPX - (P * DELTAVFEMGOPX) - (.008314 * TK * (math.log(GAMMAFEMGOPX)))) / (.008314 * TK)

            # ** math.exp ? **
        KDGAROPX = math.exp(KDGAROPX)


        if  (MODECRD > 0.01):
            #REM "CALCULATES A CORRECTED KD(GRT-CRD) "
            DELTAHFEMGCRD = (((.5 * HCRD) + ((1 / 3) * HALM)) - ((.5 * HFECRD) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGCRD = (((.5 * SCRD) + ((1 / 3) * SALM)) - ((.5 * SFECRD) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGCRD = ((.5 * VCRD) + ((1 / 3) * VALM)) - ((.5 * VFECRD) + ((1 / 3) * VPY))
            GAMMAFEMGCRD = GAMMAGAR * GAMMACRD
            KDGARCRD = ((TK * DELTASFEMGCRD) - DELTAHFEMGCRD - (P * DELTAVFEMGCRD) - (.008314 * TK * (math.log(GAMMAFEMGCRD)))) / (.008314 * TK)
            KDGARCRD = math.exp(KDGARCRD)

        #END if

        if  (MODEBT > 0.01):
            #REM "CALCULATES A CORRECTED KD(GRT-BT)"
            DELTAHFEMGBT = ((((1 / 3) * HPHL) + ((1 / 3) * HALM)) - (((1 / 3) * HANN) + ((1 / 3) * HPY))) / 1000
            DELTASFEMGBT = ((((1 / 3) * SPHL) + ((1 / 3) * SALM)) - (((1 / 3) * SANN) + ((1 / 3) * SPY))) / 1000
            DELTAVFEMGBT = (((1 / 3) * VPHL) + ((1 / 3) * VALM)) - (((1 / 3) * VANN) + ((1 / 3) * VPY))
            GAMMAGARBT = GAMMABT * GAMMAGAR
            KDGARBT = ((TK * DELTASFEMGBT) - DELTAHFEMGBT - (P * DELTAVFEMGBT) - (.008314 * TK * (math.log(GAMMAGARBT)))) / (.008314 * TK)
            KDGARBT = math.exp(KDGARBT)
        #END if

        #REM "QUADRATIC SOLUTION TO CORRECTED MG-RATIOS OF MINERALS"
        #REM "REQUIRES 10 ITERATIONS (L = 1 TO 10)."

        for L in range(10):

            #REM "QUADRATIC SOLUTION TO CORRECTED MG-RATIO OF GARNET AND OPX"
            A = MFOPX - (KDGAROPX * MFOPX)
            B = (MFOPX * KDGAROPX) + MFGAR + (XMGROCK * KDGAROPX) - XMGROCK - (KDGAROPX * MGRATIOCRD * MFCRD) - (KDGAROPX * MGRATIOBT * MFBT) + (MFCRD * MGRATIOCRD) + (MFBT * MGRATIOBT)
            C = (MGRATIOCRD * MFCRD * KDGAROPX) + (MGRATIOBT * MFBT * KDGAROPX) - (XMGROCK * KDGAROPX)
            MGRATIOOPX = (-B + (((B ** 2) - (4 * A * C)) ** .5)) / (2 * A)
            MGRATIOGAR = (XMGROCK - (MGRATIOOPX * MFOPX) - (MGRATIOCRD * MFCRD) - (MGRATIOBT * MFBT)) / MFGAR

            if  (MODEBT > 0.01):
                #REM "QUADRATIC SOLUTION TO CORRECTED MG-RATIO OF GARNET AND BIOTITE ""
                A = MFBT - (KDGARBT * MFBT)
                B = (MFBT * KDGARBT) + MFGAR + (XMGROCK * KDGARBT) - XMGROCK - (KDGARBT * MGRATIOOPX * MFOPX) - (KDGARBT * MGRATIOCRD * MFCRD) + (MFOPX * MGRATIOOPX) + (MFCRD * MGRATIOCRD)
                C = (MGRATIOOPX * MFOPX * KDGARBT) + (MGRATIOCRD * MFCRD * KDGARBT) - (XMGROCK * KDGARBT)
                MGRATIOBT = (-B + (((B ** 2) - (4 * A * C)) ** .5)) / (2 * A)
                MGRATIOGAR = (XMGROCK - (MGRATIOBT * MFBT) - (MGRATIOOPX * MFOPX) - (MGRATIOCRD * MFCRD)) / MFGAR
        #END if

            if  (MODECRD > 0.01):
                #REM "QUADRATIC SOLUTION TO CORRECTED MG-RATIO OF GARNET AND CORDIERITE ""
                A = MFCRD - (KDGARCRD * MFCRD)
                B = (MFCRD * KDGARCRD) + MFGAR + (XMGROCK * KDGARCRD) - XMGROCK - (KDGARCRD * MGRATIOOPX * MFOPX) - (KDGARCRD * MGRATIOBT * MFBT) + (MFOPX * MGRATIOOPX) + (MFBT * MGRATIOBT)
                C = (MGRATIOOPX * MFOPX * KDGARCRD) + (MGRATIOBT * MFBT * KDGARCRD) - (XMGROCK * KDGARCRD)
                MGRATIOCRD = (-B + (((B ** 2) - (4 * A * C)) ** .5)) / (2 * A)
                MGRATIOGAR = (XMGROCK - (MGRATIOCRD * MFCRD) - (MGRATIOOPX * MFOPX) - (MGRATIOBT * MFBT)) / MFGAR
        #NEXT L
        #END if

        FERATIOOPX = 1 - MGRATIOOPX
        XMGOPX = (MGRATIOOPX) * ((FE2OPX + MGOPX) / 2)
        XFEOPX = (1 - MGRATIOOPX) * ((FE2OPX + MGOPX) / 2)
        XMGGAR = (MGRATIOGAR) * ((FEGAR + MGGAR) / (FEGAR + MGGAR + CAGAR + MNGAR))
        XFEGAR = (1 - MGRATIOGAR) * ((FEGAR + MGGAR) / (FEGAR + MGGAR + CAGAR + MNGAR))
        XMGCRD = MGRATIOCRD * (XFECRD + XMGCRD)
        XFECRD = (1 - MGRATIOCRD) * (XFECRD + XMGCRD)
        XMGBT =  MGRATIOBT * (XFEBT + XMGBT)
        XFEBT =  (1 - MGRATIOBT) * (XFEBT + XMGBT)

    #NEXT I

    ##REM "CALCULATE GRT-OPX FE-MG T TO SEE if  AGREES WITH FINAL FE-AL-OPX T"
    KDGAROPX = (XFEGAR * XMGOPX) / (XMGGAR * XFEOPX)
    TGAROPX = (DELTAHFEMGOPX + (P * DELTAVFEMGOPX)) / (DELTASFEMGOPX - (.008314 * math.log(KDGAROPX)) - (.008314 * math.log(GAMMAFEMGOPX)))
    TGAROPX = TGAROPX - 273

    if  (MODECRD > 0.01):
        #REM "CALCULATE GRT-CRD FE-MG T TO SEE if  AGREES WITH FINAL FE-AL-OPX T"
        KDGARCRD = (XFEGAR * XMGCRD) / (XMGGAR * XFECRD)
        TGARCRD = (DELTAHFEMGCRD + (P * DELTAVFEMGCRD)) / (DELTASFEMGCRD - (.008314 * math.log(KDGARCRD)) - (.008314 * math.log(GAMMAFEMGCRD)))
        TCGARCRD = TGARCRD - 273
    #END if

    if  (MODEBT > 0.01):
        #REM "CALCULATE GRT-BT FE-MG T TO SEE if  AGREES WITH FINAL FE-AL-OPX T"
        KDGARBT = (XFEGAR * XMGBT) / (XMGGAR * XFEBT)
        TGARBT = (DELTAHFEMGBT + (P * DELTAVFEMGBT)) / (DELTASFEMGBT - (.008314 * math.log(KDGARBT)) - (.008314 * math.log(GAMMAGARBT)))
        TCGARBT = TGARBT - 273
    #END if

    #REM "CHECK if  RECALCULATED XMGROCK IS THE SAME AS THE INITIAL XMGROCKI"
    TEST = (MGRATIOGAR * MFGAR) + (MGRATIOOPX * MFOPX) + (MGRATIOCRD * MFCRD) + (MGRATIOBT * MFBT)

    #REM "CALCULATE Dif FERENCES (FINAL - INITIAL)"
    DifFEMGT = TGAROPX - TGAROPXI
    DifFEALT = TC - TFEALI
    DifFEMGP = P - PFEMGGTOPI
    DifFEALP = P - PFEALI
    DifGRTCRDT = TCGARCRD - TCGARCRDI
    DifGRTCRDP = P - PGARCRDI
    DifGRTBTP = P - PGARBTI
    DifGRTBTT = TCGARBT - TCGARBTI
    DifMFMR = TEST - XMGROCKI   # only used once
    DifMFMG = MGRATIOGAR - MGRATIOGARI
    DifMFMO = MGRATIOOPX - MGRATIOOPXI
    DifMFMB = MGRATIOBT  - MGRATIOBTI
    DifMFMC = MGRATIOCRD - MGRATIOCRDI

    #REM "print OUT RESULTS"
    # ** Change **
    # change to print out to a file



    print("")
    print("INITIAL AND CONVERGED P-T ESTIMATES AND MINERAL COMPOSITIONS")
    print("")

    print("             CONVERGED (FINAL)                       INITIAL                                Dif FERENCE")
    f.write("             CONVERGED (FINAL)                       INITIAL                                Dif FERENCE")


    # print #2, USING "Fe Al  GOPQ#####. C ##.## KB ";TC     ;P; USING "#####. C ##.## KB "; TFEALI  ; PFEALI; USING "####. C ##.## KB"; Dif FEALT; Dif FEALP
    #print( "Fe Al  GOPQ " + str(TC) + "C " + str(P) + "KB " + str(TFEALI) + "C " + str(PFEALI) + "KB " + str( DifFEALT) + "C " + str(DifFEALP) + " KB")
    print("Fe Al  GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TC, P, TFEALI, PFEALI, DifFEALT, DifFEALP))
    f.write("Fe Al  GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TC, P, TFEALI, PFEALI, DifFEALT, DifFEALP))


    #print , USING "GrtOpx GOPQ#####. C ##.## KB ";TGAROPX;P; USING "#####. C ##.## KB "; TGAROPXI; PFEMGGTOPI; USING "####. C ##.## KB"; Dif FEMGT; Dif FEMGP
    #print("GrtOpx GOPQ " + str(TGAROPX) + "C " + str(P) +"KB " + str(TGAROPXI) + "C " + str(PFEMGGTOPI) + "KB " + str(DifFEMGT) + "C "+ str(DifFEMGP) +" KB" )
    print("GrtOpx GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TGAROPX, P, TGAROPXI, PFEMGGTOPI, DifFEMGT, DifFEMGP))
    f.write("GrtOpx GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TGAROPX, P, TGAROPXI, PFEMGGTOPI, DifFEMGT, DifFEMGP))

    if  (MODEBT >  0.01):
        #print("GrtBt  GOPQ "+ str(TCGARBT) +"C "+ str(P) +"KB " + str(TCGARBTI) +"C "+ str(PGARBTI) +"KB " + str(DifGRTBTT) +"C " + str(DifGRTBTP) + " KB")
        print("GrtBt  GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TCGARBT, P, TCGARBTI, PGARBTI, DifGRTBTT, DifGRTBTP))
        f.write("GrtBt  GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TCGARBT, P, TCGARBTI, PGARBTI, DifGRTBTT, DifGRTBTP))

    if  (MODECRD > 0.01):
        #print , USING "GrtCrd GOPQ#####. C ##.## KB ";TCGARCRD;P; USING "#####. C ##.## KB "; TCGARCRDI;PGARCRDI; USING "####. C ##.## KB"; Dif GRTCRDT; Dif GRTCRDP
        #print("GrtCRD GOPQ " + str(TCGARCRD) + "C " + str(P) + "KB " + str(TCGARCRDI) + "C " + str(PGARCRDI) + "KB " + str(DifGRTCRDT) + "C " + str(DifGRTCRDP) + " KB")
        print("GrtCRD GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TCGARCRD, P, TCGARCRDI, PGARCRDI, DifGRTCRDT, DifGRTCRDP))
        f.write("GrtCRD GOPQ =  {:0.3f} C {:0.3f} KB \t {:0.3f} C {:0.3f} KB \t\t {:0.3f}C {:0.3f} KB".format(TCGARCRD, P, TCGARCRDI, PGARCRDI, DifGRTCRDT, DifGRTCRDP))

        #print , USING "M/FM rock         #.###     ";TEST     ; USING "        #.###     "; XMGROCKI        ; USING "       ##.###     "; Dif MFMR
        #print("M/FM rock        " + str(TEST) +"         " + str(XMGROCKI) + "        " + str(DifMFMR))
        print("M/FM rock =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(TEST, XMGROCKI, DifMFMR))
        f.write("M/FM rock =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(TEST, XMGROCKI, DifMFMR))

        #print , USING "M/FM Grt          #.###     ";MGRATIOGAR;USING "        #.###     "; MGRATIOGARI     ; USING "       ##.###     "; Dif MFMG
        #print("M/FM Grt         " + str(MGRATIOGAR) + "       " + str(MGRATIOGARI) + "         " + str(DifMFMG) )
        print("M/FM Grt =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOGAR, MGRATIOGARI, DifMFMG))
        f.write("M/FM Grt =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOGAR, MGRATIOGARI, DifMFMG))

        #print("M/FM Opx         " + str(MGRATIOOPX) + "         " + str(MGRATIOOPXI) + "      " + str(DifMFMO) )
        print("M/FM Opx =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOOPX, MGRATIOOPXI, DifMFMO))
        f.write("M/FM Opx =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOOPX, MGRATIOOPXI, DifMFMO))


    if  (MODEBT > 0.01):
        #print , USING "M/FM Bt           #.###     ";MGRATIOBT ;USING "        #.###     "; MGRATIOBTI      ; USING "       ##.###     "; Dif MFMB
        #print("M/FM Bt          " + str(MGRATIOBT) + "       " + str(MGRATIOBTI) + "          " + str(DifMFMB) )
        print("M/FM Bt =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOBT, MGRATIOBTI, DifMFMB))
        f.write("M/FM Bt =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOBT, MGRATIOBTI, DifMFMB))

    if  (MODECRD> 0.01):
        #print , USING "M/FM Crd          #.###     ";MGRATIOCRD;USING "        #.###     "; MGRATIOCRDI     ; USING "       ##.###     "; Dif MFMC
        #print("M/FM Crd         " + str(MGRATIOCRD) + "          " + str(MGRATIOCRDI) + "       " + str(DifMFMC) )
        print("M/FM Crd  =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOCRD, MGRATIOCRDI, DifMFMC))
        f.write("M/FM Crd  =  {:0.3f} \t {:0.3f} \t\t {:0.3f} ".format(MGRATIOCRD, MGRATIOCRDI, DifMFMC))

    print("")

    print ("")

    with open('export.csv', 'w', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(['Type', 'CONVERGED(FINAL)', '', '', 'INITIAL', '', '', 'Dif FERENCE'])
            filewriter.writerow(['Fe Al GOPQ', '{0:} C'.format(math.ceil(TC)), '{0:.2f} KB'.format(P), '', '{0:} C'.format(math.ceil(TFEALI)), '{0:.2f} KB'.format(PFEALI), '','{0:} C'.format(math.ceil(DifFEALT)), '{0:.2f} KB'.format(DifFEALP) ])
            filewriter.writerow(['GrtOpx GOPQ', '{0:} C'.format(math.ceil(TGAROPX)), '{0:.2f} KB'.format(P), '','{0:} C'.format(math.ceil(TGAROPXI)), '{0:.2f} KB'.format(PFEMGGTOPI),'','{0:} C'.format(math.ceil(DifFEMGT)), '{0:.2f} KB'.format( DifFEMGP) ])
            filewriter.writerow(['GrtBt  GOPQ', '{0:} C'.format(math.ceil(TCGARBT)), '{0:.2f} KB'.format(P), '','{0:} C'.format(math.ceil(TCGARBTI)), '{0:.2f} KB'.format(PGARBTI), '', '{0:} C'.format(math.ceil( DifGRTBTT)), '{0:.2f} KB'.format(  DifGRTBTP) ])
            filewriter.writerow(['GrtCRD GOPQ', '{0:} C'.format(math.ceil(TCGARCRD)), '{0:.2f} KB'.format(P),'', '{0:} C'.format(math.ceil(TCGARCRDI)), '{0:.2f} KB'.format(PGARCRDI),'','{0:} C'.format(math.ceil( DifGRTCRDT)), '{0:.2f} KB'.format(  DifGRTCRDP) ])

            filewriter.writerow(['M/FM rock', '{0:.3f}'.format(TEST), '','','{0:.3f}'.format(XMGROCKI),'', '','{0:.3f}'.format(DifMFMR)])
            filewriter.writerow(['M/FM Grt', '{0:.3f}'.format(MGRATIOGAR),'','', '{0:.3f}'.format(MGRATIOGARI),'', '','{0:.3f}'.format(DifMFMG)])
            filewriter.writerow(['M/FM Opx', '{0:.3f}'.format(MGRATIOOPX),'','', '{0:.3f}'.format(MGRATIOOPXI),'', '','{0:.3f}'.format(DifMFMO)])
            filewriter.writerow(['M/FM Bt', '{0:.3f}'.format(MGRATIOBT),'','', '{0:.3f}'.format(MGRATIOBTI),'', '','{0:.3f}'.format(DifMFMB)])
            filewriter.writerow(['M/FM Crd', '{0:.3f}'.format(MGRATIOCRD),'', '','{0:.3f}'.format(MGRATIOCRDI),'', '','{0:.3f}'.format(DifMFMC)])








    #INPUT "ANOTHER AL-IN-OPX MODEL (Y/y)"; ANSWER1$
    ANSWER1 = input("Another Model please enter Y/y: ")

    if  (ANSWER1 == "Y"):
        Model()
        return
    if  (ANSWER1 == "y"):
        Model()
        return

    print ("")
    #INPUT "ANOTHER SAMPLE (Y/y)"; ANSWER2$
    ANSWER2 = input("Another Sample please enter Y/y: ")

    if  (ANSWER2 == "Y"):
        Sample()
        return
    elif  (ANSWER2 == "y"):
        Sample()
        return
    else:
        sys.exit()









def Sample():

    global XFEOPX, XMGOPX, XALM1, ALCHOICE2, ALCHOICE
    global FEGAR, MNGAR, MGGAR, CAGAR, MODEGAR
    global SIOPX, TIOPX, ALOPX, CROPX, Fe3OPX, FE2OPX, MNOPX, MGOPX, CAOPX, MODEOPX
    global FECRD, MNCRD, MGCRD, MODECRD
    global SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT, MODEBT
    global CA, NA, KPL
    global ANSWER1, TITLE



    ANSWER1 = "Z"

    #fname = input('Enter a filename: ')
    #Title = input('Type in a Title or Short Description of the file: ')

    # $inputs are strings, else int

    #with open(fname) as f:
    #    content = f.readlines()


    # Reading in Inputs from somewhere
    # ** Need Clarification **


    print("\n PCFM-1 \n")

    filename = input("Enter file name: ")
    #readFile(filename)

    file = open(filename, "r")
    lines = file.read().splitlines()

    dataList = []
    for line in lines:
        if( lines.index(line) % 2 == 1):
            mylist = line.replace(' ','').split(',')
            mylist.pop(0)
            dataList.append(mylist)

    FEGAR, MNGAR, MGGAR, CAGAR, MODEGAR = float(dataList[0][0]), float(dataList[0][1]), float(dataList[0][2]), float(dataList[0][3]), float(dataList[0][4])
    SIOPX, TIOPX, ALOPX, CROPX, Fe3OPX, FE2OPX, MNOPX, MGOPX, CAOPX, MODEOPX = float(dataList[1][0]), float(dataList[1][1]), float(dataList[1][2]), float(dataList[1][3]), float(dataList[1][4]), float(dataList[1][5]), float(dataList[1][6]), float(dataList[1][7]), float(dataList[1][8]), float(dataList[1][9])
    FECRD, MNCRD, MGCRD, MODECRD = float(dataList[2][0]), float(dataList[2][1]), float(dataList[2][2]), float(dataList[2][3])
    SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT, MODEBT = float(dataList[3][0]), float(dataList[3][1]), float(dataList[3][2]), float(dataList[3][3]), float(dataList[3][4]), float(dataList[3][5]), float(dataList[3][6]), float(dataList[3][7]), float(dataList[3][8])
    CA, NA, KPL = float(dataList[4][0]), float(dataList[4][1]), float(dataList[4][2])

    # Print("\n\n 165-241 \n\n")
    # FEGAR = input('Enter Fe (1.842 otherwise): ')
    # MNGAR = input('Enter Mn (0.175 otherwise): ')
    # MGGAR = input('Enter Mg (0.849 otherwise): ')
    # CAGAR = input('Enter Ca (0.191 otherwise): ')
    # MODEGAR = input('Enter Mode (10.0 otherwise): ')
    #
    # SIOPX = input('Enter Si (1.864 otherwise): ')
    # TIOPX = input('Enter Ti (0.003 otherwise): ')
    # ALOPX = input('Enter Al (0.227 otherwise): ')
    # CROPX = input('Enter Cr (0.000 otherwise): ')
    # Fe3OPX = input('Enter Fe3 (0.037 otherwise): ')
    # FE2OPX = input('Enter Fe2 (0.851 otherwise): ')
    # MNOPX = input('Enter Mn (0.029 otherwise): ')
    # MGOPX = input('Enter Cr (0.969 otherwise): ')
    # CAOPX = input('Enter Ca (0.015 otherwise): ')
    # MODEOPX = input('Enter Cr (40.0 otherwise): ')
    #
    # FECRD = input('Enter Fe (0.000 otherwise): ')
    # MNCRD = input('Enter Mn (0.000 otherwise): ')
    # MGCRD = input('Enter Mg (0.000 otherwise): ')
    # MODECRD = input('Enter Mode (0.0 otherwise): ')
    #
    # SIBT = input('Enter Si (0.000 otherwise): ')
    # TIBT = input('Enter Ti (0.000 otherwise): ')
    # ALBT = input('Enter Al (0.000 otherwise): ')
    # FEBT = input('Enter Cr (0.000 otherwise): ')
    # MNBT = input('Enter Mn (0.000 otherwise): ')
    # MGBT = input('Enter Cr (0.000 otherwise): ')
    # NABT = input('Enter Ca (0.000 otherwise): ')
    # KBT = input('Enter K (0.000 otherwise): ')
    # MODEBT = input('Enter Cr (0.000 otherwise): ')
    #
    # CA = input('Enter Ca (0.354 otherwise): ')
    # NA = input('Enter Na (0.581 otherwise): ')
    # KPL = input('Enter K (0.045 otherwise): ')




    TITLE = input("Type in a title or short description of the Sample: ")

    Model()
#END of Program




#DATA -5265317.1, 341.5824, 621.4269, -3287.931, -15081040, 2211865100
#DATA -6284733.7, 268.8, 590.9042, -2826.956, -13320810, 1260328500
#DATA -6632861, 255.15, 573.43042, -2039.405, -18887168, 2319311872
#DATA -4228730, 200.1861, 439.36938, -3734.149, 0, -317023232
#DATA -908626.8, 44.2068, 80.01199, -240.276, -3546684, 491568384
#DATA -1546037.1, 66.18, 166.5795, -1200.588, -2270560, 279150300
#DATA -1192860, 96.5587, 174.2024, -1392.959, -454390, -37711400
#DATA -1631665.9, 35.375, 119.38, 774.808, -6509130, 422877600
#DATA -6216676.7, 325.9239, 610.37988, -2083.781, -21533008, 2841040896
#DATA -5155234.4, 405.01, 727.208, -4775.04, -13831900, 2119060000
#DATA -9161425.7, 416.2714, 954.3865, -7962.274, -2317258, -370214090
#DATA -8429860.2, 482.8282, 983.479, -8403.659, -1870290, -85683500

# always the same
DATASET = [[-5265317.1, 341.5824, 621.4269, -3287.931, -15081040, 2211865100], \
[-6284733.7, 268.8, 590.9042, -2826.956, -13320810, 1260328500], \
[-6632861, 255.15, 573.43042, -2039.405, -18887168, 2319311872], \
[-4228730, 200.1861, 439.36938, -3734.149, 0, -317023232], \
[-908626.8, 44.2068, 80.01199, -240.276, -3546684, 491568384], \
[-1546037.1, 66.18, 166.5795, -1200.588, -2270560, 279150300], \
[-1192860, 96.5587, 174.2024, -1392.959, -454390, -37711400], \
[ -1631665.9, 35.375, 119.38, 774.808, -6509130, 422877600], \
[-6216676.7, 325.9239, 610.37988, -2083.781, -21533008, 2841040896], \
[-5155234.4, 405.01, 727.208, -4775.04, -13831900, 2119060000], \
[-9161425.7, 416.2714, 954.3865, -7962.274, -2317258, -370214090], \
[-8429860.2, 482.8282, 983.479, -8403.659, -1870290, -85683500]]





# open DataFile
# String split on comma
# DataSet[j][i] = ?


#j, i = 12, 6
#DataSet = [[0 for x in range(12)] for y in range(6)]
#for x in range(0,j):
#	for y in range(0,i):
#	    DataSet[j][i] =

print("RCLC - GRT-OPX-PL-QTZ P-T ESTIMATES CORRECTED FOR LATE FE-MG EXCHANGE")
print("\n")


# GLOBAL VARIABLES #
FEGAR = 1.550
MNGAR = 0.200
MGGAR = 1.140
CAGAR = 0.110
MODEGAR = 10.0

SIOPX = 1.900
TIOPX = 0.010
ALOPX = 0.200
CROPX = 0.000
Fe3OPX = 0.000
FE2OPX = 0.660
MNOPX = 0.010
MGOPX = 1.230
CAOPX = 0.000
MODEOPX = 10.0

FECRD = 0.340
MNCRD = 0.010
MGCRD = 1.660
MODECRD = 10.0

SIBT = 2.760
TIBT = 0.240
ALBT = 1.300
FEBT = 0.780
MNBT = 0.010
MGBT = 1.780
NABT = 0.010
KBT = 0.970
MODEBT = 10.0

CA = 0.310
NA = 0.670
KPL = 0.020

XFEOPX = 0
XMGOPX = 0
XALM1 = 0
ALCHOICE = ""
ALCHOICE2 = ""


HALM = 0
HPY = 0
HGR = 0
HAN = 0
HBQ = 0
HEN = 0
HFS= 0
HALOPX = 0
HPHL = 0
HANN = 0
HCRD = 0
HFECRD = 0


SALM = 0
SPY = 0
SGR = 0
SAN = 0
SBQ = 0
SEN = 0
SFS = 0
SALOPX = 0
SPHL = 0
SANN = 0
SCRD = 0
SFECRD = 0



VALM = 0
VPY = 0
VGR = 0
VAN =0
VBN = 0
VEN =0
VFS = 0
VALOPX =0
VPHL= 0
VANN= 0
VCRD= 0
VFECRD= 0

AGR = 0
APY = 0
AAL = 0
GAMMAGAR =0
GAMMACRD =0
GAMMABT = 0
ANN = 0
AEN = 0
AFS = 0
AALOPX =0
GAMMAOPX =0

XMNGAR = 0
XFEGAR = 0
XCAGAR = 0
XMNGAR = 0
XFEGARI = 0
MGRATIOGA = 0
MGRATIOGARI = 0
MODXCAGAR = 0
FERATIOOPX = 0


TK = 0
PBARS = 0
P = 0
TITLE = ""



#
# print("{} {} {} {} {}".format(FEGAR, MNGAR, MGGAR, CAGAR, MODEGAR))
# print("{} {} {} {} {} {} {} {} {} {}".format(SIOPX, TIOPX, ALOPX, CROPX, Fe3OPX, FE2OPX, MNOPX, MGOPX, CAOPX, MODEOPX))
# print("{} {} {} {}".format(FECRD, MNCRD, MGCRD, MODECRD))
# print("{} {} {} {} {} {} {} {} {}".format(SIBT, TIBT, ALBT, FEBT, MNBT, MGBT, NABT, KBT, MODEBT))
# print("{} {} {}".format(CA, NA, KPL))

Sample()


# Sample:

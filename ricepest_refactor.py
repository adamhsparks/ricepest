#NAME: RICEPEST Spatial Model
#DESCRIPTION: Generic RICEPEST tool for calculating attainable and actual yields based on supplied parameters
#REQUIREMENTS: ArcGIS Spatial Analyst Extension
#DEVELOPED BY: Confidence Duku (AfricaRice), Adam Sparks (IRRI), Sander Zwart (AfricaRice
#BASED ON WORK DONE BY: Laetitia Willocquet and Serge Savary (IRRI)

#ORIGINAL VERSION of RICEPEST WITH UNMODIFIED PS3

import georasters
import numpy
import os
from src.model_funcs import (calc_actual_radiation, 
                            calc_init_temp_sum_seed, 
                            calc_init_temp_sum_transplants, 
                            calc_sum_temp_above_TBASE,
                            WEEDS,
                            WH,
                            SHR,
                            LAI,
                            BS,
                            SLA,
                            RDIST,
                            PART,
                            RDIST,
                            SHB,
                            damage_mech_bacterial_leaf_blight,
                            damage_mech_rice_leaf_blast,


                            
                            )

arcpy.AddMessage("\n            NAME:             RICEPEST Spatial Model")
arcpy.AddMessage("              DEVELOPED BY:     Confidence Duku (AfricaRice), Adam Sparks (IRRI) and Sander Zwart (AfricaRice)")
arcpy.AddMessage("              BASED ON WORK BY: Laetitia Willocquet and Serge Savary (IRRI)")
arcpy.AddMessage("              REQUIREMENTS:     ArcGIS Spatial Analyst Extension")


arcpy.AddMessage("\nSetting Environment Variables")

#obtaining directory for script and associated folders
scriptPath = sys.argv[0]
PathName = os.path.dirname(scriptPath)

#Setting environment variables;
path = PathName + "\\Output\\"
path1 = PathName + "\\Yields\\"

#analysisType = sys.argv[4]
analysisType = "Actual Yield"
cropEst = "Transplanted"
TransplantDays = 20
prodSituation = "PS2"

arcpy.AddMessage("RUNNING, " + prodSituation + " ," + analysisType)
arcpy.CheckOutExtension("Spatial")

#Prefix the climate and disease data with the following to distinguish them
RadGDB = arcpy.ListRasters("rad*", "GRID")
TempGDB = arcpy.ListRasters("tmean*", "GRID")
IniTemp = arcpy.ListRasters("i*", "GRID")
BlastGDB = arcpy.ListRasters("*blast*", "TIF")
BlightGDB = arcpy.ListRasters("*bblight*", "TIF")



# -------------------------------
if prodSituation == "PS1":
    PANW = 0
    LEAFW = 10
    STEMW = 6
    ROOTW = 5
    TBASE = 8

    #calculating the sum of temperature between crop establishment and the start of simulation
    jt = 0
    rt = 0
    kt = 0

    #check if rice variety is direct seeded or transplanted.
    if cropEst == "Direct Seeded":
        arcpy.AddMessage("\nCalculating Sum of Initial Temperature for Direct Seeded Systems")
        rt = calc_init_temp_sum_seed(path, IniTemp, TBASE, rt)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        out2 = calc_init_temp_sum_transplants(path, TransplantDays, IniTemp, TBASE, jt, rt, kt)
    else:
        raise ValueError("Invalid crop establishment method: " + cropEst, expects = "Direct Seeded or Transplanted")
    
    #calculating sum of temp above TBASE
    calc_sum_temp_above_TBASE(TempGDB, TBASE, out2)

    #------------------ Damage Mechanisms ------------------
    #! Whats this????? 1 x 0 ??????
        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        
        blbrf = damage_mech_bacterial_leaf_blight(analysisType, oneRaster, BlightGDB, tempDate )

            #Damage mechanism due to Rice Leaf Blast
        rlbrf = damage_mech_rice_leaf_blast(analysisType, oneRaster, BlastGDB, tempDate, blast_start_date_idx=6  )

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        
        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

     
    this_RAD = calc_actual_radiation(RadGDB, tempDate)


 
    # ----------- Calculating coefficients of partitioning ---------------
    arcpy.AddMessage("  Calculating coefficients of partitioning")

        #Increase in dry weight of organs
        this_rdist = RDIST()
        this_part = PART()
        this_partp = this_part[1]
        this_partl = this_part[0]
        this_partst = this_part[2]
        this_partr = this_part[3]
        this_rrsenl = RRSENL()
        this_shb = SHB()
        this_trshbl = this_shb[0]

        arcpy.AddMessage("  Calculating the dry weight of organs\n")

        #Calculate the increase in dry weight for panicles
        ppanw1 = arcpy.sa.Plus(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = arcpy.sa.Times(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = arcpy.sa.Times(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = arcpy.sa.Minus(this_partst, this_rdist)

        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = arcpy.sa.Times(this_rrsenl, LEAFW)
        else:
            rsenl = arcpy.sa.Times(this_rrsenl, thisLEAFW)

        pLEAFW = arcpy.sa.Minus(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = arcpy.sa.Minus(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = arcpy.sa.Plus(ppanw3, PANW)
            thisSTEMW = arcpy.sa.Plus(pstemw, STEMW)
            thisROOTW = arcpy.sa.Plus(ROOTW, this_partr)
            thisLEAFW = arcpy.sa.Plus(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)

    output =  path1 + "Yield_PS1"
    thisPANW.save(output)
    arcpy.AddMessage("COMPLETED")


elif prodSituation == "PS2":
    PANW = 0
    LEAFW = 6
    STEMW = 4
    ROOTW = 3
    TBASE = 8

    #calculating the sum of temperature between crop establishment and the start of simulation
    jt = 0
    rt = 0
    kt = 0

    #check if rice variety is directed seeded or transplanted.
    if cropEst == "Direct Seeded":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Direct Seeded Systems")
        rt = calc_init_temp_sum_seed(path, IniTemp, TBASE, rt)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        out2 = calc_init_temp_sum_transplants(path, TransplantDays, IniTemp, TBASE, jt, rt, kt)
    else:
        raise ValueError("Invalid crop establishment method: " + cropEst, expects = "Direct Seeded or Transplanted")
    


    #calculating sum of temp above TBASE

    calc_sum_temp_above_TBASE(TempGDB, TBASE, out2)

        
#! Whats this????? 1 x 0 ??
        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        
        blbrf = damage_mech_bacterial_leaf(analysisType, oneRaster, BlightGDB, tempDate )

    
        rlbrf = damage_mech_rice_leaf_blast(analysisType, oneRaster, BlastGDB, tempDate, blast_start_date_idx=6  )
        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)



        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        

        

this_RAD = calc_actual_radiation(RadGDB, tempDate)


        #Calculating the Rate of Growth and the Pool of assimilates
        def POOL (k=-0.6):
            arcpy.AddMessage("  Calculating Pool of assimilates")
            where1 = "VALUE <= 0.9"
            RUE = arcpy.sa.Con(co_sumt, 1.1, 1, where1)
            this_LAI = LAI()
            pRG1 = arcpy.sa.Times(k, this_LAI)
            pRG2 = arcpy.sa.Exp(pRG1)
            pRG3 = arcpy.sa.Minus(1, pRG2)
            pRG4 = arcpy.sa.Times(pRG3, this_RAD)
            pRG5 = arcpy.sa.Times(pRG4, RUE)

            #Calculating the rate of growth after treatment with Weeds
            this_weedrf = WEEDS()
            RG = arcpy.sa.Times(pRG5, this_weedrf)

            #mRG = arcpy.sa.Times(RG, Days)
            ty = str(i)
            name = "pool" + ty
            RG.save(path + name)
            return RG


        arcpy.AddMessage("  Calculating coefficients of partitioning")
        #Calculating the coefficient of partitioning of leaves relative to DVS
        

        #Calculating the coefficient of partitioning of roots relative to DVS
        

        #Calculating the coefficient of partitioning of panicles relative to DVS
        


        

        

        

        #Calculating the relative rate of senescence
        def RRSENL ():
            arcpy.AddMessage("  Calculating the relative rate of senescence")
            #sumt107 = (-425 * (2.07**2))+(2550 * 2.07)-2125
            sumt13 = (-435 * (2.3**2))+(2755 * 2.3)-2320
            #sumt163 = (-425 * (2.63**2))+(2550 * 2.63)-2125
            where1 = "\"VALUE\" <= %d" % sumt13
            #where2 = "\"VALUE\" <= %d" % sumt13
            #where3 = "\"VALUE\" <= %d" % sumt163
            rrsenl = arcpy.sa.Con(co_sumt, 0, 0.04, where1)
            return rrsenl

        #Increase in dry weight of organs
        this_rdist = RDIST()
        this_part = PART()
        this_partp = this_part[1]
        this_partl = this_part[0]
        this_partst = this_part[2]
        this_partr = this_part[3]
        this_rrsenl = RRSENL()
        this_shb = SHB()
        this_trshbl = this_shb[0]

        arcpy.AddMessage("  Calculating the dry weight of organs\n")

        #Calculate the increase in dry weight for panicles
        ppanw1 = arcpy.sa.Plus(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = arcpy.sa.Times(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = arcpy.sa.Times(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = arcpy.sa.Minus(this_partst, this_rdist)


        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = arcpy.sa.Times(this_rrsenl, LEAFW)
        else:
            rsenl = arcpy.sa.Times(this_rrsenl, thisLEAFW)

        pLEAFW = arcpy.sa.Minus(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = arcpy.sa.Minus(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = arcpy.sa.Plus(ppanw3, PANW)
            thisSTEMW = arcpy.sa.Plus(pstemw, STEMW)
            thisROOTW = arcpy.sa.Plus(ROOTW, this_partr)
            thisLEAFW = arcpy.sa.Plus(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)

    output =  path1 + "Yield_PS2"
    thisPANW.save(output)
    arcpy.AddMessage("COMPLETED")


elif prodSituation == "PS3":
    PANW = 0
    LEAFW = 17
    STEMW = 15
    ROOTW = 7
    TBASE = 8

    #calculating the sum of temperature between crop establishment and the start of simulation
    jt = 0
    rt = 0
    kt = 0

    #check if rice variety is direct seeded or transplanted.
    if cropEst == "Direct Seeded":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Direct Seeded Systems")
        rt = calc_init_temp_sum_seed(path, IniTemp, TBASE, rt)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        out2 = calc_init_temp_sum_transplants(path, TransplantDays, IniTemp, TBASE, jt, rt, kt)
    else:
        raise ValueError("Invalid crop establishment method: " + cropEst, expects = "Direct Seeded or Transplanted")
    


    #calculating sum of temp above TBASE

    calc_sum_temp_above_TBASE(TempGDB, TBASE, out2)

        def sla():
            arcpy.AddMessage("  Calculating Specific Leaf Area")
            #calculating SLA
            #always add 1 to the value of sumt**
            #hence sumt007 with value of 0.07 becomes 1.07 in the equation
            sumt0 = (-650 * (1**2))+(3750 * 1)-3100
            sumt1 = (-650 * (2**2))+(3750 * 2)-3100
            sumt2 = (-650 * (3**2))+(3750 * 3)-3100
            where0 = "\"VALUE\" <= %d" % sumt0
            where1 = "\"VALUE\" <= %d" % sumt1
            where2 = "\"VALUE\" <= %d" % sumt2
            SLA = arcpy.sa.Con(co_sumt, 0.037, (arcpy.sa.Con(co_sumt, 0.018, 0.017, where1)), where0)
            arcpy.BuildPyramids_management(SLA)
            return SLA
#! Whats this????? 1 x 0 ??
        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)


##        #Damage mechanism due to Brown Spot
##        if analysisType == "Attainable Yield":
##            brownSpot = 0
##            bsbrf = arcpy.sa.Times(oneRaster, brownSpot)
##
##        elif analysisType == "Actual Yield":
##            arcpy.AddMessage("  Calculating Reduction Factor for BS")
##            for brownSpot in BSpotGDB:
##                bsDate = brownSpot[2:]
##                if tempDate == bsDate:
##                    pbsrf = arcpy.sa.Divide(bsdm, 100)
##                    pbsrf1 = arcpy.sa.Minus(1, pbsbrf)
##                    bsbrf = arcpy.sa.Power(pbsbrf1, 6.3)
##                    break
##            else:
##                bsbrf = 1

        blbrf = damage_mech_bacterial_leaf(analysisType, oneRaster, BlightGDB, tempDate )

        rlbrf = damage_mech_rice_leaf_blast(analysisType, oneRaster, BlastGDB, tempDate, blast_start_date_idx=6  )


        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)



        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        

        

this_RAD = calc_actual_radiation(RadGDB, tempDate)


        #Calculating the Rate of Growth and the Pool of assimilates
        def POOL (k=-0.6):
            arcpy.AddMessage("  Calculating Pool of assimilates")
            sumt09 = (-650 * (1.9**2))+(3750 * 1.9)-3100
            where09 = "\"VALUE\" <= %d" % sumt09
            RUE = arcpy.sa.Con(co_sumt, 1.4, 0.8, where09)
            this_LAI = LAI()
            pRG1 = arcpy.sa.Times(k, this_LAI)
            pRG2 = arcpy.sa.Exp(pRG1)
            pRG3 = arcpy.sa.Minus(1, pRG2)
            pRG4 = arcpy.sa.Times(pRG3, this_RAD)
            pRG5 = arcpy.sa.Times(pRG4, RUE)

            #Calculating the rate of growth after treatment with Weeds
            this_weedrf = WEEDS()
            RG = arcpy.sa.Times(pRG5, this_weedrf)
            #mRG = arcpy.sa.Times(RG, Days)
            ty = str(i)
            name = "pool" + ty
            RG.save(path + name)
            return RG


        arcpy.AddMessage("  Calculating coefficients of partitioning")
        #Calculating the coefficient of partitioning of leaves relative to DVS
        
        #! WHERE DO THESE GET USED?????
        def COPARTL ():
            sumt00 = (-650 * (1**2))+(3750 * 1)-3100
            sumt05 = (-650 * (1.5**2))+(3750 * 1.5)-3100
            where0 = "\"VALUE\" <= %d" % sumt00
            where05 = "\"VALUE\" <= %d" % sumt05
            copartl = arcpy.sa.Con(co_sumt, 0.55, (arcpy.sa.Con(co_sumt, 0.45, 0, where05)), where0)
            return copartl

        #Calculating the coefficient of partitioning of roots relative to DVS
        def COPARTR ():
            sumt00 = (-650 * (1**2))+(3750 * 1)-3100
            sumt08 = (-650 * (1.8**2))+(3750 * 1.8)-3100
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            copartr = arcpy.sa.Con(co_sumt, 0.3, arcpy.sa.Con(co_sumt, 0.1, 0, where08), where0)
            return copartr

        #Calculating the coefficient of partitioning of panicles relative to DVS
        def COPARTP ():
            sumt075 = (-650 * (1.75**2))+(3750 * 1.75)-3100
            sumt11 = (-650 * (2.1**2))+(3750 * 2.1)-3100
            where075 = "\"VALUE\" > %d" % sumt075
            where11 = "\"VALUE\" >= %d" % sumt11
            copartp = arcpy.sa.Con(co_sumt, 1, arcpy.sa.Con(co_sumt, 0.3, 0, where075), where11)
            return copartp




        

        #Calculating the relative rate of senescence
        def RRSENL ():
            arcpy.AddMessage("  Calculating the relative rate of senescence")
            #sumt107 = (-425 * (2.07**2))+(2550 * 2.07)-2125
            sumt13 = (-650 * (2.3**2))+(3750 * 2.3)-3100
            #sumt163 = (-425 * (2.63**2))+(2550 * 2.63)-2125
            where13 = "\"VALUE\" >= %d" % sumt13
            #where2 = "\"VALUE\" <= %d" % sumt13
            #where3 = "\"VALUE\" <= %d" % sumt163
            rrsenl = arcpy.sa.Con(co_sumt, 0.05, 0, where13)
            return rrsenl

        #Increase in dry weight of organs
        this_rdist = RDIST(co_sumt, param1 = -650, param2 = 3750, param3 = 3100)
        this_part = PART()
        this_partp = this_part[1]
        this_partl = this_part[0]
        this_partst = this_part[2]
        this_partr = this_part[3]
        this_rrsenl = RRSENL()
        this_shb = SHB()
        this_trshbl = this_shb[0]

        arcpy.AddMessage("  Calculating the dry weight of organs\n")
        #Calculate the increase in dry weight for panicles
        ppanw1 = arcpy.sa.Plus(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = arcpy.sa.Times(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = arcpy.sa.Times(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = arcpy.sa.Minus(this_partst, this_rdist)

        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = arcpy.sa.Times(this_rrsenl, LEAFW)
        else:
            rsenl = arcpy.sa.Times(this_rrsenl, thisLEAFW)

        pLEAFW = arcpy.sa.Minus(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = arcpy.sa.Minus(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = arcpy.sa.Plus(ppanw3, PANW)
            thisSTEMW = arcpy.sa.Plus(pstemw, STEMW)
            thisROOTW = arcpy.sa.Plus(ROOTW, this_partr)
            thisLEAFW = arcpy.sa.Plus(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)

    output =  path1 + "yield_ps3"
    thisPANW.save(output)
    arcpy.AddMessage("COMPLETED")


elif prodSituation == "PS4":

    PANW = 0
    LEAFW = 15
    STEMW = 15
    ROOTW = 10
    TBASE = 8

    #calculating the sum of temperature between crop establishment and the start of simulation
    jt = 0
    rt = 0
    kt = 0

    #check if rice variety is directed seeded or transplanted.
    if cropEst == "Direct Seeded":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Direct Seeded Systems")
        rt = calc_init_temp_sum_seed(path, IniTemp, TBASE, rt)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        out2 = calc_init_temp_sum_transplants(path, TransplantDays, IniTemp, TBASE, jt, rt, kt)
        #! lines below this: out22 = 0.785 * out4 indented within else statement in PS4.
    else:
        raise ValueError("Invalid crop establishment method: " + cropEst, expects = "Direct Seeded or Transplanted")
    
    #calculating sum of temp above TBASE
    calc_sum_temp_above_TBASE(TempGDB, TBASE, out2)

        def sla():
            arcpy.AddMessage("  Calculating Specific Leaf Area")
            #calculating SLA
            #always add 1 to the value of sumt**
            #hence sumt007 with value of 0.07 becomes 1.07 in the equation
            sumt0 = (-300 * (1**2))+(2000 * 1)-1700
            sumt05 = (-300 * (1.5**2))+(2000 * 1.5)-1700
            sumt1 = (-300 * (2**2))+(2000 * 2)-1700
            sumt2 = (-300 * (3**2))+(2000 * 3)-1700
            where0 = "\"VALUE\" <= %d" % sumt0
            where05 = "\"VALUE\" <= %d" % sumt05
            where1 = "\"VALUE\" <= %d" % sumt1
            where2 = "\"VALUE\" <= %d" % sumt2
            SLA = arcpy.sa.Con(co_sumt, 0.035, (arcpy.sa.Con(co_sumt, 0.027, 0.022, where05)), where0)
            arcpy.BuildPyramids_management(SLA)
            return SLA
#! Whats this????? 1 x 0 ??
        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)
        
        #Damage mechanism due to Bacterial Leaf Blight
        if analysisType == "Attainable Yield":
            blight = 0
            blbrf = arcpy.sa.Minus(oneRaster, blight)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLIGHT")
            for blight in BlightGDB:
                blightDate = blight[25:]
                if tempDate == blightDate:
                    pblbrf = arcpy.sa.Divide(blight, 100)
                    blbrf = arcpy.sa.Minus(1, pblbrf)
                    break
            else:
                blbrf = 1

        #Damage mechanism due to Rice Leaf Blast
        #! idx changes to 25 ?????
        rlbrf = damage_mech_rice_leaf_blast(analysisType, oneRaster, BlastGDB, tempDate, blast_start_date_idx=25  )

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)



        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        

        

this_RAD = calc_actual_radiation(RadGDB, tempDate)

        #Calculating the Rate of Growth and the Pool of assimilates
        def POOL (k=-0.6):
            arcpy.AddMessage("  Calculating Pool of assimilates")
            where1 = "VALUE <= 0.9"
            RUE = arcpy.sa.Con(co_sumt, 1.3, 1.2, where1)
            this_LAI = LAI()
            pRG1 = arcpy.sa.Times(k, this_LAI)
            pRG2 = arcpy.sa.Exp(pRG1)
            pRG3 = arcpy.sa.Minus(1, pRG2)
            pRG4 = arcpy.sa.Times(pRG3, this_RAD)
            pRG5 = arcpy.sa.Times(pRG4, RUE)

            #Calculating the rate of growth after treatment with Weeds
            this_weedrf = WEEDS()
            RG = arcpy.sa.Times(pRG5, this_weedrf)
            #mRG = arcpy.sa.Times(RG, Days)
            ty = str(i)
            name = "pool" + ty
            RG.save(path + name)
            return RG

        arcpy.AddMessage("  Calculating coefficients of partitioning")
        #Calculating the coefficient of partitioning of leaves relative to DVS
        def COPARTL ():
            sumt00 = (-300 * (1**2))+(2000 * 1)-1700
            sumt05 = (-300 * (1.5**2))+(2000 * 1.5)-1700
            where0 = "\"VALUE\" <= %d" % sumt00
            where05 = "\"VALUE\" <= %d" % sumt05
            copartl = arcpy.sa.Con(co_sumt, 0.55, (arcpy.sa.Con(co_sumt, 0.45, 0, where05)), where0)
            return copartl

        #Calculating the coefficient of partitioning of roots relative to DVS
        def COPARTR ():
            sumt00 = (-300 * (1**2))+(2000 * 1)-1700
            sumt08 = (-300 * (1.8**2))+(2000 * 1.8)-1700
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            copartr = arcpy.sa.Con(co_sumt, 0.3, arcpy.sa.Con(co_sumt, 0.1, 0, where08), where0)
            return copartr

        #Calculating the coefficient of partitioning of panicles relative to DVS
        def COPARTP ():
            sumt065 = (-300 * (1.65**2))+(2000 * 1.65)-1700
            sumt11 = (-300 * (2.1**2))+(2000 * 2.1)-1700
            where065 = "\"VALUE\" > %d" % sumt065
            where11 = "\"VALUE\" >= %d" % sumt11
            copartp = arcpy.sa.Con(co_sumt, 1, arcpy.sa.Con(co_sumt, 0.3, 0, where065), where11)
            return copartp




        

        



        #Calculating the relative rate of senescence
        def RRSENL ():
            arcpy.AddMessage("  Calculating the relative rate of senescence")
            #sumt107 = (-425 * (2.07**2))+(2550 * 2.07)-2125
            sumt13 = (-300 * (2.3**2))+(2000 * 2.3)-1700
            #sumt163 = (-425 * (2.63**2))+(2550 * 2.63)-2125
            where1 = "\"VALUE\" <= %d" % sumt13
            #where2 = "\"VALUE\" <= %d" % sumt13
            #where3 = "\"VALUE\" <= %d" % sumt163
            rrsenl = arcpy.sa.Con(co_sumt, 0, 0.04, where1)
            return rrsenl

        #Increase in dry weight of organs
        this_rdist = RDIST(co_sumt, param1 = -300, param2 = 2000, param3 = 1700)
        this_part = PART()
        this_partp = this_part[1]
        this_partl = this_part[0]
        this_partst = this_part[2]
        this_partr = this_part[3]
        this_rrsenl = RRSENL()
        this_shb = SHB()
        this_trshbl = this_shb[0]

        arcpy.AddMessage("  Calculating the dry weight of organs\n")

        #Calculate the increase in dry weight for panicles
        ppanw1 = arcpy.sa.Plus(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = arcpy.sa.Times(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = arcpy.sa.Times(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = arcpy.sa.Minus(this_partst, this_rdist)

        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = arcpy.sa.Times(this_rrsenl, LEAFW)
        else:
            rsenl = arcpy.sa.Times(this_rrsenl, thisLEAFW)

        pLEAFW = arcpy.sa.Minus(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = arcpy.sa.Minus(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = arcpy.sa.Plus(ppanw3, PANW)
            thisSTEMW = arcpy.sa.Plus(pstemw, STEMW)
            thisROOTW = arcpy.sa.Plus(ROOTW, this_partr)
            thisLEAFW = arcpy.sa.Plus(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)

    output =  path1 + "Yield_PS4"
    thisPANW.save(output)
    arcpy.AddMessage("COMPLETED")

elif prodSituation == "PS5":
    PANW = 0
    LEAFW = 8
    STEMW = 5
    ROOTW = 4.5
    TBASE = 8

    #calculating the sum of temperature between crop establishment and the start of simulation
    jt = 0
    rt = 0
    kt = 0

    #check if rice variety is directed seeded or transplanted.
    if cropEst == "Direct Seeded":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Direct Seeded Systems")
        rt = calc_init_temp_sum_seed(path, IniTemp, TBASE, rt)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        out2 = calc_init_temp_sum_transplants(path, TransplantDays, IniTemp, TBASE, jt, rt, kt)
    else:
        raise ValueError("Invalid crop establishment method: " + cropEst, expects = "Direct Seeded or Transplanted")
    


    #calculating sum of temp above TBASE

    calc_sum_temp_above_TBASE(TempGDB, TBASE, out2)

        
#! Whats this????? 1 x 0 ??
        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)
        
        blbrf = damage_mech_bacterial_leaf(analysisType, oneRaster, BlightGDB, tempDate )

        #Damage mechanism due to Rice Leaf Blast
        rlbrf = damage_mech_rice_leaf_blast(analysisType, oneRaster, BlastGDB, tempDate, blast_start_date_idx=6  )

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)



        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        

        

this_RAD = calc_actual_radiation(RadGDB, tempDate)


        #Calculating the Rate of Growth and the Pool of assimilates
        def POOL (k=-0.6):
            arcpy.AddMessage("  Calculating Pool of assimilates")
            where1 = "VALUE <= 0.9"
            RUE = arcpy.sa.Con(co_sumt, 1.2, 1.1, where1)
            this_LAI = LAI()
            pRG1 = arcpy.sa.Times(k, this_LAI)
            pRG2 = arcpy.sa.Exp(pRG1)
            pRG3 = arcpy.sa.Minus(1, pRG2)
            pRG4 = arcpy.sa.Times(pRG3, this_RAD)
            pRG5 = arcpy.sa.Times(pRG4, RUE)

            #Calculating the rate of growth after treatment with Weeds
            this_weedrf = WEEDS()
            RG = arcpy.sa.Times(pRG5, this_weedrf)
            #mRG = arcpy.sa.Times(RG, Days)
            ty = str(i)
            name = "pool" + ty
            RG.save(path + name)
            return RG


        arcpy.AddMessage("  Calculating coefficients of partitioning")
        #Calculating the coefficient of partitioning of leaves relative to DVS
        

        #Calculating the coefficient of partitioning of roots relative to DVS
        

        #Calculating the coefficient of partitioning of panicles relative to DVS
        

        #Calculating the relative rate of senescence
        def RRSENL ():
            arcpy.AddMessage("  Calculating the relative rate of senescence")
            #sumt107 = (-425 * (2.07**2))+(2550 * 2.07)-2125
            sumt13 = (-435 * (2.3**2))+(2755 * 2.3)-2320
            #sumt163 = (-425 * (2.63**2))+(2550 * 2.63)-2125
            where1 = "\"VALUE\" <= %d" % sumt13
            #where2 = "\"VALUE\" <= %d" % sumt13
            #where3 = "\"VALUE\" <= %d" % sumt163
            rrsenl = arcpy.sa.Con(co_sumt, 0, 0.06, where1)
            return rrsenl

        #Increase in dry weight of organs
        this_rdist = RDIST()
        this_part = PART()
        this_partp = this_part[1]
        this_partl = this_part[0]
        this_partst = this_part[2]
        this_partr = this_part[3]
        this_rrsenl = RRSENL()
        this_shb = SHB()
        this_trshbl = this_shb[0]

        arcpy.AddMessage("  Calculating the dry weight of organs\n")

        #Calculate the increase in dry weight for panicles
        ppanw1 = arcpy.sa.Plus(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = arcpy.sa.Times(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = arcpy.sa.Times(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = arcpy.sa.Minus(this_partst, this_rdist)

        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = arcpy.sa.Times(this_rrsenl, LEAFW)
        else:
            rsenl = arcpy.sa.Times(this_rrsenl, thisLEAFW)

        pLEAFW = arcpy.sa.Minus(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = arcpy.sa.Minus(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = arcpy.sa.Plus(ppanw3, PANW)
            thisSTEMW = arcpy.sa.Plus(pstemw, STEMW)
            thisROOTW = arcpy.sa.Plus(ROOTW, this_partr)
            thisLEAFW = arcpy.sa.Plus(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)

    output =  path1 + "Yield_PS5"
    thisPANW.save(output)
    arcpy.AddMessage("COMPLETED")


elif prodSituation == "PS6":
    PANW = 0
    LEAFW = 8
    STEMW = 5
    ROOTW = 4.5
    TBASE = 8

    #calculating the sum of temperature between crop establishment and the start of simulation
    jt = 0
    rt = 0
    kt = 0

    #check if rice variety is direct seeded or transplanted.
    if cropEst == "Direct Seeded":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Direct Seeded Systems")
        rt = calc_init_temp_sum_seed(path, IniTemp, TBASE, rt)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        out2 = calc_init_temp_sum_transplants(path, TransplantDays, IniTemp, TBASE, jt, rt, kt)
    else:
        raise ValueError("Invalid crop establishment method: " + cropEst, expects = "Direct Seeded or Transplanted")
    


    #calculating sum of temp above TBASE

    calc_sum_temp_above_TBASE(TempGDB, TBASE, out2)

        
#! Whats this????? 1 x 0 ??
        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        
        blbrf = damage_mech_bacterial_leaf(analysisType, oneRaster, BlightGDB, tempDate )

        #Damage mechanism due to Rice Leaf Blast
        rlbrf = damage_mech_rice_leaf_blast(analysisType, oneRaster, BlastGDB, tempDate, blast_start_date_idx=6  )


        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)



        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        

        

this_RAD = calc_actual_radiation(RadGDB, tempDate)


        #Calculating the Rate of Growth and the Pool of assimilates
        def POOL (k=-0.6):
            arcpy.AddMessage("  Calculating Pool of assimilates")
            where1 = "VALUE <= 0.9"
            RUE = arcpy.sa.Con(co_sumt, 1.3, 1.2, where1)
            this_LAI = LAI()
            pRG1 = arcpy.sa.Times(k, this_LAI)
            pRG2 = arcpy.sa.Exp(pRG1)
            pRG3 = arcpy.sa.Minus(1, pRG2)
            pRG4 = arcpy.sa.Times(pRG3, this_RAD)
            pRG5 = arcpy.sa.Times(pRG4, RUE)

            #Calculating the rate of growth after treatment with Weeds
            this_weedrf = WEEDS()
            RG = arcpy.sa.Times(pRG5, this_weedrf)
            #mRG = arcpy.sa.Times(RG, Days)
            ty = str(i)
            name = "pool" + ty
            RG.save(path + name)
            return RG


        arcpy.AddMessage("  Calculating coefficients of partitioning")
        #Calculating the coefficient of partitioning of leaves relative to DVS
        

        #Calculating the coefficient of partitioning of roots relative to DVS
        

        #Calculating the coefficient of partitioning of panicles relative to DVS
        



        

        

        

        #Calculating the relative rate of senescence
        def RRSENL ():
            arcpy.AddMessage("  Calculating the relative rate of senescence")
            #sumt107 = (-425 * (2.07**2))+(2550 * 2.07)-2125
            sumt13 = (-435 * (2.3**2))+(2755 * 2.3)-2320
            #sumt163 = (-425 * (2.63**2))+(2550 * 2.63)-2125
            where1 = "\"VALUE\" <= %d" % sumt13
            #where2 = "\"VALUE\" <= %d" % sumt13
            #where3 = "\"VALUE\" <= %d" % sumt163
            rrsenl = arcpy.sa.Con(co_sumt, 0, 0.06, where1)
            return rrsenl

        #Increase in dry weight of organs
        this_rdist = RDIST()

        this_part = PART()
        this_partp = this_part[1]
        this_partl = this_part[0]
        this_partst = this_part[2]
        this_partr = this_part[3]
        this_rrsenl = RRSENL()
        this_shb = SHB()
        this_trshbl = this_shb[0]

        arcpy.AddMessage("  Calculating the dry weight of organs\n")
        #Calculate the increase in dry weight for panicles
        ppanw1 = arcpy.sa.Plus(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = arcpy.sa.Times(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = arcpy.sa.Times(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = arcpy.sa.Minus(this_partst, this_rdist)


        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = arcpy.sa.Times(this_rrsenl, LEAFW)
        else:
            rsenl = arcpy.sa.Times(this_rrsenl, thisLEAFW)

        pLEAFW = arcpy.sa.Minus(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = arcpy.sa.Minus(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = arcpy.sa.Plus(ppanw3, PANW)
            thisSTEMW = arcpy.sa.Plus(pstemw, STEMW)
            thisROOTW = arcpy.sa.Plus(ROOTW, this_partr)
            thisLEAFW = arcpy.sa.Plus(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(path + name)

    output =  path1 + "Yield_PS6"
    thisPANW.save(output)
    arcpy.AddMessage("COMPLETED")

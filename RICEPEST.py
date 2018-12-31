#NAME: RICEPEST Spatial Model
#DESCRIPTION: Generic RICEPEST tool for calculating attainable and actual yields based on supplied parameters
#REQUIREMENTS: ArcGIS Spatial Analyst Extension
#DEVELOPED BY: Confidence Duku (AfricaRice), Adam Sparks (IRRI), Sander Zwart (AfricaRice
#BASED ON WORK DONE BY: Laetitia Willocquet and Serge Savary (IRRI)

#ORIGINAL VERSION WITH UNMODIFIED PS3

import arcpy
import sys
import os

arcpy.AddMessage("\n            NAME:             RICEPEST Spatial Model")
arcpy.AddMessage("              DEVELOPED BY:     Confidence Duku (AfricaRice), Adam Sparks (IRRI) and Sander Zwart (AfricaRice)")
arcpy.AddMessage("              BASED ON WORK BY: Laetitia Willocquet and Serge Savary (IRRI)")
arcpy.AddMessage("              REQUIREMENTS:     ArcGIS Spatial Analyst Extension")


arcpy.AddMessage("\nSetting Environment Variables")

#obtaining directory for script and associated folders
scriptPath = sys.argv[0]
PathName = os.path.dirname(scriptPath)

#Setting environment variables;
arcpy.env.extent = "MINOF"
arcpy.env.scratchWorkspace = PathName + "\\Output"
arcpy.env.cellSize = "MINOF"
arcpy.env.workspace =  PathName + "\\ab50"
arcpy.env.overwriteOutput = 1
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
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out2 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out2 = out2 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if jt == 0:
                if kt != TransplantDays:
                    out3 = arcpy.sa.Minus(raster, TBASE)
                    out4 = arcpy.sa.Con(out3, 0, out3, where1)
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break
            else:
                if kt != TransplantDays:
                    out1 = arcpy.sa.Minus(raster, TBASE)
                    out6 = arcpy.sa.Con(out1, 0, out1, where1)
                    out4 = out4 + out6
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break

        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out23 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                out23.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out23 = out23 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                name1 = "sumtr"
                out23.save(path + name)

        out22 = 0.785 * out4
        out22.save(path + "tshock")

        out2 = out23 - out22
        out2.save(path + name1)


    #calculating sum of temp above TBASE

    i = 0
    for temp in TempGDB:
        test = 0
        where1 = "\"VALUE\" < %d" % test
        tempDate = temp[5:]
        #calculating sum of temperature
        bn = i + 15
        ft = str(bn)
        arcpy.AddMessage("SIMULATION AT " + ft + " DACE")
        arcpy.AddMessage("  Calculating Sum of Temperature above TBASE")
        if i == 0:
            out4 = arcpy.sa.Minus(temp, TBASE)
            out5 = arcpy.sa.Con(out4, 0, out4, where1)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)
        else:
            DTemp = arcpy.sa.Minus(temp, TBASE)
            DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)

        def sla():
            arcpy.AddMessage("  Calculating Specific Leaf Area")
            #calculating SLA
            #always add 1 to the value of sumt**
            #hence sumt007 with value of 0.07 becomes 1.07 in the equation
            sumt0 = (-435 * (1**2))+(2755 * 1)-2320
            sumt05 = (-435 * (1.5**2))+(2755 * 1.5)-2320
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            sumt2 = (-435 * (3**2))+(2755 * 3)-2320
            where0 = "\"VALUE\" <= %d" % sumt0
            where05 = "\"VALUE\" <= %d" % sumt05
            where1 = "\"VALUE\" <= %d" % sumt1
            where2 = "\"VALUE\" <= %d" % sumt2
            SLA = arcpy.sa.Con(co_sumt, 0.037, (arcpy.sa.Con(co_sumt, 0.03, 0.017, where05)), where0)
            arcpy.BuildPyramids_management(SLA)
            return SLA

        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        def SHB ():
            #The first Damage mechanism due to Sheath Blight
            prshbl = arcpy.sa.Times(0.00076, Sevl)
            if i == 1:
                rshbl = arcpy.sa.Times(prshbl, LEAFW)
            else:
                rshbl = arcpy.sa.Times(prshbl, thisLEAFW)
            #trshbl = arcpy.sa.Times(rshbl, Days)

            #The second Damage mechanism due to Sheath Blight
            pshbrf = arcpy.sa.Divide(Sevl, 100)
            shbrf = arcpy.sa.Minus(1, pshbrf)
            return [rshbl, shbrf]


        #Damage mechanism due to Brown Spot
        mybsdm = 0
        bsdm = arcpy.sa.Times(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = arcpy.sa.Divide(bsdm, 100)
            pbsrf1 = arcpy.sa.Minus(1, pbsrf)
            bsrf = arcpy.sa.Power(pbsrf1, beta)
            return bsrf

        #Damage mechanism due to Bacterial Leaf Blight
        if analysisType == "Attainable Yield":
            blight = 0
            blbrf = arcpy.sa.Minus(oneRaster, blight)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLIGHT")
            for blight in BlightGDB:
                blightDate = blight[7:]
                if tempDate == blightDate:
                    pblbrf = arcpy.sa.Divide(blight, 100)
                    blbrf = arcpy.sa.Minus(1, pblbrf)
                    break
            else:
                blbrf = 1

        #Damage mechanism due to Rice Leaf Blast
        if analysisType == "Attainable Yield":
            blast = 0
            rlbrf = arcpy.sa.Minus(oneRaster, blast)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLAST")
            for blast in BlastGDB:
                blastDate = blast[6:]
                if tempDate == blastDate:
                    prlbrf = arcpy.sa.Divide(blast, 100)
                    prlbrf2 = arcpy.sa.Minus(1, prlbrf)
                    rlbrf = arcpy.sa.Power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        def SHR ():
            shrrf = arcpy.sa.Minus(1, shrdm)
            #shrrf = arcpy.sa.Times(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        def WH ():
            whrf = arcpy.sa.Minus(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = arcpy.sa.Minus(1, prfwd1)
            weedrf = arcpy.sa.Minus(1, rfwd)
            return weedrf

        #Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight
        def LAI ():
            arcpy.AddMessage("  Calculating Leaf Area Index")
            thisSLA = sla()
            this_shb1 = SHB()
            this_shbrf = this_shb1[1]
            this_bsrf = BS()

            if i == 1:
                pLAI = arcpy.sa.Times(thisSLA, LEAFW)
            else:
                pLAI = arcpy.sa.Times(thisSLA, thisLEAFW)

            pLAI1 = arcpy.sa.Times(blbrf, pLAI)
            pLAI2 = arcpy.sa.Times(rlbrf, pLAI1)
            pLAI3 = arcpy.sa.Times(this_shbrf, pLAI2)
            LAI = arcpy.sa.Times(this_bsrf, pLAI3)
            return LAI

        #Calculating actual radiation
        for rad in RadGDB:
            radDate = rad[3:]
            if tempDate == radDate:
                this_RAD = rad
                break


        #Calculating the Rate of Growth and the Pool of assimilates
        def POOL (k=-0.6):
            arcpy.AddMessage("  Calculating Pool of assimilates")
            sumt09 = (-435 * (1.9**2))+(2755 * 1.9)-2320
            where09 = "VALUE <= 0.9"
            RUE = arcpy.sa.Con(co_sumt, 1.3, 1.2, where09)
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
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt07 = (-435 * (1.7**2))+(2755 * 1.7)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where07 = "\"VALUE\" <= %d" % sumt07
            copartl = arcpy.sa.Con(co_sumt, 0.55, (arcpy.sa.Con(co_sumt, 0.45, 0, where07)), where0)
            return copartl

        #Calculating the coefficient of partitioning of roots relative to DVS
        def COPARTR ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt08 = (-435 * (1.8**2))+(2755 * 1.8)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            copartr = arcpy.sa.Con(co_sumt, 0.3, arcpy.sa.Con(co_sumt, 0.1, 0, where08), where0)
            return copartr

        #Calculating the coefficient of partitioning of panicles relative to DVS
        def COPARTP ():
            sumt075 = (-435 * (1.75**2))+(2755 * 1.75)-2320
            sumt11 = (-435 * (2.1**2))+(2755 * 2.1)-2320
            where075 = "\"VALUE\" > %d" % sumt075
            where11 = "\"VALUE\" >= %d" % sumt11
            copartp = arcpy.sa.Con(co_sumt, 1, arcpy.sa.Con(co_sumt, 0.3, 0, where075), where11)
            return copartp

        #Calculating the coefficient of partitioning of stems relative to DVS
        this_copartp = COPARTP ()
        this_copartl = COPARTL ()
        def COPARTST ():
            pr_copartst = arcpy.sa.Minus(this_copartl, this_copartp)
            copartst = arcpy.sa.Minus(1, pr_copartst)
            return copartst

        #Partitioning of assimilates
        this_copartst = COPARTST()
        this_POOL = POOL()

        def PART ():
            arcpy.AddMessage("  Calculating partitioning of assimilates")
            #Calculate the partitioning of assimilates to leaves
            this_copartr = COPARTR()
            this = arcpy.sa.Minus(1, this_copartr)
            pr_partl = arcpy.sa.Times(this_POOL, this_copartl)
            partl = arcpy.sa.Times(pr_partl, this)

            #Calculate the partitioning of assimilates to panicles
            pr_partp = arcpy.sa.Times(this_POOL, this_copartp)
            partp = arcpy.sa.Times(pr_partp, this)

            #Calculate the partitioning of assimilates to stems
            pr_partst = arcpy.sa.Times(this_POOL, this_copartst)
            partst = arcpy.sa.Times(pr_partst, this)

            #Calculate the partitioning of assimilates to roots
            partr = arcpy.sa.Times(this_POOL, this_copartr)

            return[partl, partp, partst, partr]

        #Redistribution of reserves accumulated in the stems
        def RDIST():
            arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            where1 = "\"VALUE\" < %d" % sumt1
            rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
            return rdist

        #Calculating the relative rate of senescence
        def RRSENL ():
            arcpy.AddMessage("  Calculating the relative rate of senescence")
            #sumt107 = (-425 * (2.07**2))+(2550 * 2.07)-2125
            sumt13 = (-435 * (2.3**2))+(2755 * 2.3)-2320
            #sumt163 = (-425 * (2.63**2))+(2550 * 2.63)-2125
            where1 = "\"VALUE\" < %d" % sumt13
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
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out2 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out2 = out2 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if jt == 0:
                if kt != TransplantDays:
                    out3 = arcpy.sa.Minus(raster, TBASE)
                    out4 = arcpy.sa.Con(out3, 0, out3, where1)
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break
            else:
                if kt != TransplantDays:
                    out1 = arcpy.sa.Minus(raster, TBASE)
                    out6 = arcpy.sa.Con(out1, 0, out1, where1)
                    out4 = out4 + out6
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break

        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out23 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                out23.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out23 = out23 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                name1 = "sumtr"
                out23.save(path + name)

        out22 = 0.785 * out4
        out22.save(path + "tshock")

        out2 = out23 - out22
        out2.save(path + name1)


    #calculating sum of temp above TBASE

    i = 0
    for temp in TempGDB:
        test = 0
        where1 = "\"VALUE\" < %d" % test
        tempDate = temp[5:]
        #calculating sum of temperature
        bn = i + 15
        ft = str(bn)
        arcpy.AddMessage("SIMULATION AT " + ft + " DACE")
        arcpy.AddMessage("  Calculating Sum of Temperature above TBASE")
        if i == 0:
            out4 = arcpy.sa.Minus(temp, TBASE)
            out5 = arcpy.sa.Con(out4, 0, out4, where1)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)
        else:
            DTemp = arcpy.sa.Minus(temp, TBASE)
            DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)

        def sla():
            arcpy.AddMessage("  Calculating Specific Leaf Area")
            #calculating SLA
            #always add 1 to the value of sumt**
            #hence sumt007 with value of 0.07 becomes 1.07 in the equation
            sumt0 = (-435 * (1**2))+(2755 * 1)-2320
            sumt05 = (-435 * (1.5**2))+(2755 * 1.5)-2320
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            sumt2 = (-435 * (3**2))+(2755 * 3)-2320
            where0 = "\"VALUE\" <= %d" % sumt0
            where05 = "\"VALUE\" <= %d" % sumt05
            where1 = "\"VALUE\" <= %d" % sumt1
            where2 = "\"VALUE\" <= %d" % sumt2
            SLA = arcpy.sa.Con(co_sumt, 0.037, (arcpy.sa.Con(co_sumt, 0.03, 0.017, where05)), where0)
            arcpy.BuildPyramids_management(SLA)
            return SLA

        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        def SHB ():
            #The first Damage mechanism due to Sheath Blight
            prshbl = arcpy.sa.Times(0.00076, Sevl)
            if i == 1:
                rshbl = arcpy.sa.Times(prshbl, LEAFW)
            else:
                rshbl = arcpy.sa.Times(prshbl, thisLEAFW)
            #trshbl = arcpy.sa.Times(rshbl, Days)

            #The second Damage mechanism due to Sheath Blight
            pshbrf = arcpy.sa.Divide(Sevl, 100)
            shbrf = arcpy.sa.Minus(1, pshbrf)
            return [rshbl, shbrf]


        #Damage mechanism due to Brown Spot
        mybsdm = 0
        bsdm = arcpy.sa.Times(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = arcpy.sa.Divide(bsdm, 100)
            pbsrf1 = arcpy.sa.Minus(1, pbsrf)
            bsrf = arcpy.sa.Power(pbsrf1, beta)
            return bsrf

        #Damage mechanism due to Bacterial Leaf Blight
        if analysisType == "Attainable Yield":
            blight = 0
            blbrf = arcpy.sa.Minus(oneRaster, blight)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLIGHT")
            for blight in BlightGDB:
                blightDate = blight[7:]
                if tempDate == blightDate:
                    pblbrf = arcpy.sa.Divide(blight, 100)
                    blbrf = arcpy.sa.Minus(1, pblbrf)
                    break
            else:
                blbrf = 1

        #Damage mechanism due to Rice Leaf Blast
        if analysisType == "Attainable Yield":
            blast = 0
            rlbrf = arcpy.sa.Minus(oneRaster, blast)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLAST")
            for blast in BlastGDB:
                blastDate = blast[6:]
                if tempDate == blastDate:
                    prlbrf = arcpy.sa.Divide(blast, 100)
                    prlbrf2 = arcpy.sa.Minus(1, prlbrf)
                    rlbrf = arcpy.sa.Power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        def SHR ():
            shrrf = arcpy.sa.Minus(1, shrdm)
            #shrrf = arcpy.sa.Times(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        def WH ():
            whrf = arcpy.sa.Minus(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = arcpy.sa.Minus(1, prfwd1)
            weedrf = arcpy.sa.Minus(1, rfwd)
            return weedrf

        #Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight
        def LAI ():
            arcpy.AddMessage("  Calculating Leaf Area Index")
            thisSLA = sla()
            this_shb1 = SHB()
            this_shbrf = this_shb1[1]
            this_bsrf = BS()

            if i == 1:
                pLAI = arcpy.sa.Times(thisSLA, LEAFW)
            else:
                pLAI = arcpy.sa.Times(thisSLA, thisLEAFW)

            pLAI1 = arcpy.sa.Times(blbrf, pLAI)
            pLAI2 = arcpy.sa.Times(rlbrf, pLAI1)
            pLAI3 = arcpy.sa.Times(this_shbrf, pLAI2)
            LAI = arcpy.sa.Times(this_bsrf, pLAI3)
            return LAI

        #Calculating actual radiation
        for rad in RadGDB:
            radDate = rad[3:]
            if tempDate == radDate:
                this_RAD = rad
                break


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
        def COPARTL ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt07 = (-435 * (1.7**2))+(2755 * 1.7)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where07 = "\"VALUE\" <= %d" % sumt07
            copartl = arcpy.sa.Con(co_sumt, 0.55, (arcpy.sa.Con(co_sumt, 0.45, 0, where07)), where0)
            return copartl

        #Calculating the coefficient of partitioning of roots relative to DVS
        def COPARTR ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt08 = (-435 * (1.8**2))+(2755 * 1.8)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            copartr = arcpy.sa.Con(co_sumt, 0.3, arcpy.sa.Con(co_sumt, 0.1, 0, where08), where0)
            return copartr

        #Calculating the coefficient of partitioning of panicles relative to DVS
        def COPARTP ():
            sumt075 = (-435 * (1.75**2))+(2755 * 1.75)-2320
            sumt11 = (-435 * (2.1**2))+(2755 * 2.1)-2320
            where075 = "\"VALUE\" > %d" % sumt075
            where11 = "\"VALUE\" >= %d" % sumt11
            copartp = arcpy.sa.Con(co_sumt, 1, arcpy.sa.Con(co_sumt, 0.3, 0, where075), where11)
            return copartp
        #Calculating the coefficient of partitioning of stems relative to DVS
        this_copartp = COPARTP ()
        this_copartl = COPARTL ()
        def COPARTST ():
            pr_copartst = arcpy.sa.Minus(this_copartl, this_copartp)
            copartst = arcpy.sa.Minus(1, pr_copartst)
            return copartst

        #Partitioning of assimilates
        this_copartst = COPARTST()
        this_POOL = POOL()

        def PART ():
            arcpy.AddMessage("  Calculating partitioning of assimilates")
            #Calculate the partitioning of assimilates to leaves
            this_copartr = COPARTR()
            this = arcpy.sa.Minus(1, this_copartr)
            pr_partl = arcpy.sa.Times(this_POOL, this_copartl)
            partl = arcpy.sa.Times(pr_partl, this)

            #Calculate the partitioning of assimilates to panicles
            pr_partp = arcpy.sa.Times(this_POOL, this_copartp)
            partp = arcpy.sa.Times(pr_partp, this)

            #Calculate the partitioning of assimilates to stems
            pr_partst = arcpy.sa.Times(this_POOL, this_copartst)
            partst = arcpy.sa.Times(pr_partst, this)

            #Calculate the partitioning of assimilates to roots
            partr = arcpy.sa.Times(this_POOL, this_copartr)

            return[partl, partp, partst, partr]

        #Redistribution of reserves accumulated in the stems
        def RDIST():
            arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            where1 = "\"VALUE\" < %d" % sumt1
            rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
            return rdist

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
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out2 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out2 = out2 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if jt == 0:
                if kt != TransplantDays:
                    out3 = arcpy.sa.Minus(raster, TBASE)
                    out4 = arcpy.sa.Con(out3, 0, out3, where1)
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break
            else:
                if kt != TransplantDays:
                    out1 = arcpy.sa.Minus(raster, TBASE)
                    out6 = arcpy.sa.Con(out1, 0, out1, where1)
                    out4 = out4 + out6
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break

        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out23 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                out23.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out23 = out23 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                name1 = "sumtr"
                out23.save(path + name)

        out22 = 0.785 * out4
        out22.save(path + "tshock")

        out2 = out23 - out22
        out2.save(path + name1)


    #calculating sum of temp above TBASE

    i = 0
    for temp in TempGDB:
        test = 0
        where1 = "\"VALUE\" < %d" % test
        tempDate = temp[5:]
        #calculating sum of temperature
        bn = i + 15
        ft = str(bn)
        arcpy.AddMessage("SIMULATION AT " + ft + " DACE")
        arcpy.AddMessage("  Calculating Sum of Temperature above TBASE")
        if i == 0:
            out4 = arcpy.sa.Minus(temp, TBASE)
            out5 = arcpy.sa.Con(out4, 0, out4, where1)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)
        else:
            DTemp = arcpy.sa.Minus(temp, TBASE)
            DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)

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

        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        def SHB ():
            #The first Damage mechanism due to Sheath Blight
            prshbl = arcpy.sa.Times(0.00076, Sevl)
            if i == 1:
                rshbl = arcpy.sa.Times(prshbl, LEAFW)
            else:
                rshbl = arcpy.sa.Times(prshbl, thisLEAFW)
            #trshbl = arcpy.sa.Times(rshbl, Days)

            #The second Damage mechanism due to Sheath Blight
            pshbrf = arcpy.sa.Divide(Sevl, 100)
            shbrf = arcpy.sa.Minus(1, pshbrf)
            return [rshbl, shbrf]

        mybsdm = 0
        bsdm = arcpy.sa.Times(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = arcpy.sa.Divide(bsdm, 100)
            pbsrf1 = arcpy.sa.Minus(1, pbsrf)
            bsrf = arcpy.sa.Power(pbsrf1, beta)
            return bsrf

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

        #Damage mechanism due to Bacterial Leaf Blight
        if analysisType == "Attainable Yield":
            blight = 0
            blbrf = arcpy.sa.Minus(oneRaster, blight)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLIGHT")
            for blight in BlightGDB:
                blightDate = blight[7:]
                if tempDate == blightDate:
                    pblbrf = arcpy.sa.Divide(blight, 100)
                    blbrf = arcpy.sa.Minus(1, pblbrf)
                    break
            else:
                blbrf = 1

        #Damage mechanism due to Rice Leaf Blast
        if analysisType == "Attainable Yield":
            blast = 0
            rlbrf = arcpy.sa.Minus(oneRaster, blast)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLAST")
            for blast in BlastGDB:
                blastDate = blast[6:]
                if tempDate == blastDate:
                    prlbrf = arcpy.sa.Divide(blast, 100)
                    prlbrf2 = arcpy.sa.Minus(1, prlbrf)
                    rlbrf = arcpy.sa.Power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1


        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        def SHR ():
            shrrf = arcpy.sa.Minus(1, shrdm)
            #shrrf = arcpy.sa.Times(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        def WH ():
            whrf = arcpy.sa.Minus(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = arcpy.sa.Minus(1, prfwd1)
            weedrf = arcpy.sa.Minus(1, rfwd)
            return weedrf

        #Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight
        def LAI ():
            arcpy.AddMessage("  Calculating Leaf Area Index")
            thisSLA = sla()
            this_shb1 = SHB()
            this_shbrf = this_shb1[1]
            this_bsrf = BS()

            if i == 1:
                pLAI = arcpy.sa.Times(thisSLA, LEAFW)
            else:
                pLAI = arcpy.sa.Times(thisSLA, thisLEAFW)

            pLAI1 = arcpy.sa.Times(blbrf, pLAI)
            pLAI2 = arcpy.sa.Times(rlbrf, pLAI1)
            pLAI3 = arcpy.sa.Times(this_shbrf, pLAI2)
            LAI = arcpy.sa.Times(this_bsrf, pLAI3)
            return LAI

        #Calculating actual radiation
        for rad in RadGDB:
            radDate = rad[3:]
            if tempDate == radDate:
                this_RAD = rad
                break


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


        #Calculating the coefficient of partitioning of stems relative to DVS
        this_copartp = COPARTP ()
        this_copartl = COPARTL ()
        def COPARTST ():
            pr_copartst = arcpy.sa.Minus(this_copartl, this_copartp)
            copartst = arcpy.sa.Minus(1, pr_copartst)
            return copartst

        #Partitioning of assimilates
        this_copartst = COPARTST()
        this_POOL = POOL()

        def PART ():
            arcpy.AddMessage("  Calculating partitioning of assimilates")
            #Calculate the partitioning of assimilates to leaves
            this_copartr = COPARTR()
            this = arcpy.sa.Minus(1, this_copartr)
            pr_partl = arcpy.sa.Times(this_POOL, this_copartl)
            partl = arcpy.sa.Times(pr_partl, this)

            #Calculate the partitioning of assimilates to panicles
            pr_partp = arcpy.sa.Times(this_POOL, this_copartp)
            partp = arcpy.sa.Times(pr_partp, this)

            #Calculate the partitioning of assimilates to stems
            pr_partst = arcpy.sa.Times(this_POOL, this_copartst)
            partst = arcpy.sa.Times(pr_partst, this)

            #Calculate the partitioning of assimilates to roots
            partr = arcpy.sa.Times(this_POOL, this_copartr)

            return[partl, partp, partst, partr]

        #Redistribution of reserves accumulated in the stems
        def RDIST():
            arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-650 * (2**2))+(3750 * 2)-3100
            where1 = "\"VALUE\" < %d" % sumt1
            rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
            return rdist

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
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            raster = arcpy.sa.Times(raster, 1)
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, 8)
                out2 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out2 = out2 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if jt == 0:
                if kt != TransplantDays:
                    out3 = arcpy.sa.Minus(raster, TBASE)
                    out4 = arcpy.sa.Con(out3, 0, out3, where1)
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break
            else:
                if kt != TransplantDays:
                    out1 = arcpy.sa.Minus(raster, TBASE)
                    out6 = arcpy.sa.Con(out1, 0, out1, where1)
                    out4 = out4 + out6
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break

        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out23 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                out23.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out23 = out23 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                name1 = "sumtr"
                out23.save(path + name)

                out22 = 0.785 * out4
                out22.save(path + "tshock")

                out2 = out23 - out22
                out2.save(path + name1)


    #calculating sum of temp above TBASE

    i = 0
    for temp in TempGDB:
        test = 0
        where1 = "\"VALUE\" < %d" % test
        tempDate = temp[5:]
        #calculating sum of temperature
        bn = i + 15
        ft = str(bn)
        arcpy.AddMessage("SIMULATION AT " + ft + " DACE")
        arcpy.AddMessage("  Calculating Sum of Temperature above TBASE")
        if i == 0:
            out4 = arcpy.sa.Minus(temp, TBASE)
            out5 = arcpy.sa.Con(out4, 0, out4, where1)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)
        else:
            DTemp = arcpy.sa.Minus(temp, TBASE)
            DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)

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

        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        def SHB ():
            #The first Damage mechanism due to Sheath Blight
            prshbl = arcpy.sa.Times(0.00076, Sevl)
            if i == 1:
                rshbl = arcpy.sa.Times(prshbl, LEAFW)
            else:
                rshbl = arcpy.sa.Times(prshbl, thisLEAFW)
            #trshbl = arcpy.sa.Times(rshbl, Days)

            #The second Damage mechanism due to Sheath Blight
            pshbrf = arcpy.sa.Divide(Sevl, 100)
            shbrf = arcpy.sa.Minus(1, pshbrf)
            return [rshbl, shbrf]

        #Damage mechanism due to Brown Spot
        mybsdm = 0
        bsdm = arcpy.sa.Times(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = arcpy.sa.Divide(bsdm, 100)
            pbsrf1 = arcpy.sa.Minus(1, pbsrf)
            bsrf = arcpy.sa.Power(pbsrf1, beta)
            return bsrf

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
        if analysisType == "Attainable Yield":
            blast = 0
            rlbrf = arcpy.sa.Minus(oneRaster, blast)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLAST")
            for blast in BlastGDB:
                blastDate = blast[25:]
                if tempDate == blastDate:
                    prlbrf = arcpy.sa.Divide(blast, 100)
                    prlbrf2 = arcpy.sa.Minus(1, prlbrf)
                    rlbrf = arcpy.sa.Power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        def SHR ():
            shrrf = arcpy.sa.Minus(1, shrdm)
            #shrrf = arcpy.sa.Times(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        def WH ():
            whrf = arcpy.sa.Minus(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = arcpy.sa.Minus(1, prfwd1)
            weedrf = arcpy.sa.Minus(1, rfwd)
            return weedrf

        #Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight
        def LAI ():
            arcpy.AddMessage("  Calculating Leaf Area Index")
            thisSLA = sla()
            this_shb1 = SHB()
            this_shbrf = this_shb1[1]
            this_bsrf = BS()

            if i == 1:
                pLAI = arcpy.sa.Times(thisSLA, LEAFW)
            else:
                pLAI = arcpy.sa.Times(thisSLA, thisLEAFW)

            pLAI1 = arcpy.sa.Times(blbrf, pLAI)
            pLAI2 = arcpy.sa.Times(rlbrf, pLAI1)
            pLAI3 = arcpy.sa.Times(this_shbrf, pLAI2)
            LAI = arcpy.sa.Times(this_bsrf, pLAI3)
            return LAI

        #Calculating actual radiation
        for rad in RadGDB:
            radDate = rad[3:]
            if tempDate == radDate:
                this_RAD = rad
                break

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


        #Calculating the coefficient of partitioning of stems relative to DVS
        this_copartp = COPARTP ()
        this_copartl = COPARTL ()
        def COPARTST ():
            pr_copartst = arcpy.sa.Minus(this_copartl, this_copartp)
            copartst = arcpy.sa.Minus(1, pr_copartst)
            return copartst

        #Partitioning of assimilates
        this_copartst = COPARTST()
        this_POOL = POOL()

        def PART ():
            arcpy.AddMessage("  Calculating partitioning of assimilates")
            #Calculate the partitioning of assimilates to leaves
            this_copartr = COPARTR()
            this = arcpy.sa.Minus(1, this_copartr)
            pr_partl = arcpy.sa.Times(this_POOL, this_copartl)
            partl = arcpy.sa.Times(pr_partl, this)

            #Calculate the partitioning of assimilates to panicles
            pr_partp = arcpy.sa.Times(this_POOL, this_copartp)
            partp = arcpy.sa.Times(pr_partp, this)

            #Calculate the partitioning of assimilates to stems
            pr_partst = arcpy.sa.Times(this_POOL, this_copartst)
            partst = arcpy.sa.Times(pr_partst, this)

            #Calculate the partitioning of assimilates to roots
            partr = arcpy.sa.Times(this_POOL, this_copartr)

            return[partl, partp, partst, partr]

        #Redistribution of reserves accumulated in the stems
        def RDIST():
            arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-300 * (2**2))+(2000 * 2)-1700
            where1 = "\"VALUE\" < %d" % sumt1
            rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
            return rdist

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
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out2 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out2 = out2 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if jt == 0:
                if kt != TransplantDays:
                    out3 = arcpy.sa.Minus(raster, TBASE)
                    out4 = arcpy.sa.Con(out3, 0, out3, where1)
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break
            else:
                if kt != TransplantDays:
                    out1 = arcpy.sa.Minus(raster, TBASE)
                    out6 = arcpy.sa.Con(out1, 0, out1, where1)
                    out4 = out4 + out6
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break

        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out23 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                out23.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out23 = out23 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                name1 = "sumtr"
                out23.save(path + name)

        out22 = 0.785 * out4
        out22.save(path + "tshock")

        out2 = out23 - out22
        out2.save(path + name1)


    #calculating sum of temp above TBASE

    i = 0
    for temp in TempGDB:
        test = 0
        where1 = "\"VALUE\" < %d" % test
        tempDate = temp[5:]
        #calculating sum of temperature
        bn = i + 15
        ft = str(bn)
        arcpy.AddMessage("SIMULATION AT " + ft + " DACE")
        arcpy.AddMessage("  Calculating Sum of Temperature above TBASE")
        if i == 0:
            out4 = arcpy.sa.Minus(temp, TBASE)
            out5 = arcpy.sa.Con(out4, 0, out4, where1)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)
        else:
            DTemp = arcpy.sa.Minus(temp, TBASE)
            DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)

        def sla():
            arcpy.AddMessage("  Calculating Specific Leaf Area")
            #calculating SLA
            #always add 1 to the value of sumt**
            #hence sumt007 with value of 0.07 becomes 1.07 in the equation
            sumt0 = (-435 * (1**2))+(2755 * 1)-2320
            sumt05 = (-435 * (1.5**2))+(2755 * 1.5)-2320
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            sumt2 = (-435 * (3**2))+(2755 * 3)-2320
            where0 = "\"VALUE\" <= %d" % sumt0
            where05 = "\"VALUE\" <= %d" % sumt05
            where1 = "\"VALUE\" <= %d" % sumt1
            where2 = "\"VALUE\" <= %d" % sumt2
            SLA = arcpy.sa.Con(co_sumt, 0.037, (arcpy.sa.Con(co_sumt, 0.03, 0.017, where05)), where0)
            arcpy.BuildPyramids_management(SLA)
            return SLA

        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        def SHB ():
            #The first Damage mechanism due to Sheath Blight
            prshbl = arcpy.sa.Times(0.00076, Sevl)
            if i == 1:
                rshbl = arcpy.sa.Times(prshbl, LEAFW)
            else:
                rshbl = arcpy.sa.Times(prshbl, thisLEAFW)
            #trshbl = arcpy.sa.Times(rshbl, Days)

            #The second Damage mechanism due to Sheath Blight
            pshbrf = arcpy.sa.Divide(Sevl, 100)
            shbrf = arcpy.sa.Minus(1, pshbrf)
            return [rshbl, shbrf]

        #Damage mechanism due to Brown Spot
        mybsdm = 0
        bsdm = arcpy.sa.Times(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = arcpy.sa.Divide(bsdm, 100)
            pbsrf1 = arcpy.sa.Minus(1, pbsrf)
            bsrf = arcpy.sa.Power(pbsrf1, beta)
            return bsrf

        #Damage mechanism due to Bacterial Leaf Blight
        if analysisType == "Attainable Yield":
            blight = 0
            blbrf = arcpy.sa.Minus(oneRaster, blight)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLIGHT")
            for blight in BlightGDB:
                blightDate = blight[7:]
                if tempDate == blightDate:
                    pblbrf = arcpy.sa.Divide(blight, 100)
                    blbrf = arcpy.sa.Minus(1, pblbrf)
                    break
            else:
                blbrf = 1

        #Damage mechanism due to Rice Leaf Blast
        if analysisType == "Attainable Yield":
            blast = 0
            rlbrf = arcpy.sa.Minus(oneRaster, blast)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLAST")
            for blast in BlastGDB:
                blastDate = blast[6:]
                if tempDate == blastDate:
                    prlbrf = arcpy.sa.Divide(blast, 100)
                    prlbrf2 = arcpy.sa.Minus(1, prlbrf)
                    rlbrf = arcpy.sa.Power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        def SHR ():
            shrrf = arcpy.sa.Minus(1, shrdm)
            #shrrf = arcpy.sa.Times(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        def WH ():
            whrf = arcpy.sa.Minus(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = arcpy.sa.Minus(1, prfwd1)
            weedrf = arcpy.sa.Minus(1, rfwd)
            return weedrf

        #Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight
        def LAI ():
            arcpy.AddMessage("  Calculating Leaf Area Index")
            thisSLA = sla()
            this_shb1 = SHB()
            this_shbrf = this_shb1[1]
            this_bsrf = BS()

            if i == 1:
                pLAI = arcpy.sa.Times(thisSLA, LEAFW)
            else:
                pLAI = arcpy.sa.Times(thisSLA, thisLEAFW)

            pLAI1 = arcpy.sa.Times(blbrf, pLAI)
            pLAI2 = arcpy.sa.Times(rlbrf, pLAI1)
            pLAI3 = arcpy.sa.Times(this_shbrf, pLAI2)
            LAI = arcpy.sa.Times(this_bsrf, pLAI3)
            return LAI

        #Calculating actual radiation
        for rad in RadGDB:
            radDate = rad[3:]
            if tempDate == radDate:
                this_RAD = rad
                break


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
        def COPARTL ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt07 = (-435 * (1.7**2))+(2755 * 1.7)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where07 = "\"VALUE\" <= %d" % sumt07
            copartl = arcpy.sa.Con(co_sumt, 0.55, (arcpy.sa.Con(co_sumt, 0.45, 0, where07)), where0)
            return copartl

        #Calculating the coefficient of partitioning of roots relative to DVS
        def COPARTR ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt08 = (-435 * (1.8**2))+(2755 * 1.8)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            copartr = arcpy.sa.Con(co_sumt, 0.3, arcpy.sa.Con(co_sumt, 0.1, 0, where08), where0)
            return copartr

        #Calculating the coefficient of partitioning of panicles relative to DVS
        def COPARTP ():
            sumt075 = (-435 * (1.75**2))+(2755 * 1.75)-2320
            sumt11 = (-435 * (2.1**2))+(2755 * 2.1)-2320
            where075 = "\"VALUE\" > %d" % sumt075
            where11 = "\"VALUE\" >= %d" % sumt11
            copartp = arcpy.sa.Con(co_sumt, 1, arcpy.sa.Con(co_sumt, 0.3, 0, where075), where11)
            return copartp

        #Calculating the coefficient of partitioning of stems relative to DVS
        this_copartp = COPARTP ()
        this_copartl = COPARTL ()
        def COPARTST ():
            pr_copartst = arcpy.sa.Minus(this_copartl, this_copartp)
            copartst = arcpy.sa.Minus(1, pr_copartst)
            return copartst

        #Partitioning of assimilates
        this_copartst = COPARTST()
        this_POOL = POOL()

        def PART ():
            arcpy.AddMessage("  Calculating partitioning of assimilates")
            #Calculate the partitioning of assimilates to leaves
            this_copartr = COPARTR()
            this = arcpy.sa.Minus(1, this_copartr)
            pr_partl = arcpy.sa.Times(this_POOL, this_copartl)
            partl = arcpy.sa.Times(pr_partl, this)

            #Calculate the partitioning of assimilates to panicles
            pr_partp = arcpy.sa.Times(this_POOL, this_copartp)
            partp = arcpy.sa.Times(pr_partp, this)

            #Calculate the partitioning of assimilates to stems
            pr_partst = arcpy.sa.Times(this_POOL, this_copartst)
            partst = arcpy.sa.Times(pr_partst, this)

            #Calculate the partitioning of assimilates to roots
            partr = arcpy.sa.Times(this_POOL, this_copartr)

            return[partl, partp, partst, partr]

        #Redistribution of reserves accumulated in the stems
        def RDIST():
            arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            where1 = "\"VALUE\" < %d" % sumt1
            rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
            return rdist

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
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out2 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out2 = out2 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out2_" + ty
                out2.save(path + name)
    elif cropEst == "Transplanted":
        arcpy.AddMessage("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if jt == 0:
                if kt != TransplantDays:
                    out3 = arcpy.sa.Minus(raster, TBASE)
                    out4 = arcpy.sa.Con(out3, 0, out3, where1)
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break
            else:
                if kt != TransplantDays:
                    out1 = arcpy.sa.Minus(raster, TBASE)
                    out6 = arcpy.sa.Con(out1, 0, out1, where1)
                    out4 = out4 + out6
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break

        for raster in IniTemp:
            test = 0
            where1 = "\"VALUE\" < %d" % test
            if rt == 0:
                out3 = arcpy.sa.Minus(raster, TBASE)
                out23 = arcpy.sa.Con(out3, 0, out3, where1)
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                out23.save(path + name)
            else:
                out1 = arcpy.sa.Minus(raster, TBASE)
                out6 = arcpy.sa.Con(out1, 0, out1, where1)
                out23 = out23 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                name1 = "sumtr"
                out23.save(path + name)

        out22 = 0.785 * out4
        out22.save(path + "tshock")

        out2 = out23 - out22
        out2.save(path + name1)


    #calculating sum of temp above TBASE

    i = 0
    for temp in TempGDB:
        test = 0
        where1 = "\"VALUE\" < %d" % test
        tempDate = temp[5:]
        #calculating sum of temperature
        bn = i + 15
        ft = str(bn)
        arcpy.AddMessage("SIMULATION AT " + ft + " DACE")
        arcpy.AddMessage("  Calculating Sum of Temperature above TBASE")
        if i == 0:
            out4 = arcpy.sa.Minus(temp, TBASE)
            out5 = arcpy.sa.Con(out4, 0, out4, where1)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)
        else:
            DTemp = arcpy.sa.Minus(temp, TBASE)
            DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(path + name)

        def sla():
            arcpy.AddMessage("  Calculating Specific Leaf Area")
            #calculating SLA
            #always add 1 to the value of sumt**
            #hence sumt007 with value of 0.07 becomes 1.07 in the equation
            sumt0 = (-435 * (1**2))+(2755 * 1)-2320
            sumt05 = (-435 * (1.5**2))+(2755 * 1.5)-2320
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            sumt2 = (-435 * (3**2))+(2755 * 3)-2320
            where0 = "\"VALUE\" <= %d" % sumt0
            where05 = "\"VALUE\" <= %d" % sumt05
            where1 = "\"VALUE\" <= %d" % sumt1
            where2 = "\"VALUE\" <= %d" % sumt2
            SLA = arcpy.sa.Con(co_sumt, 0.037, (arcpy.sa.Con(co_sumt, 0.03, 0.017, where05)), where0)
            arcpy.BuildPyramids_management(SLA)
            return SLA

        mysevl = 0
        trysevl = co_sumt
        oneRaster = arcpy.sa.Divide(trysevl, trysevl)
        Sevl = arcpy.sa.Times(oneRaster, mysevl)

        def SHB ():
            #The first Damage mechanism due to Sheath Blight
            prshbl = arcpy.sa.Times(0.00076, Sevl)
            if i == 1:
                rshbl = arcpy.sa.Times(prshbl, LEAFW)
            else:
                rshbl = arcpy.sa.Times(prshbl, thisLEAFW)
            #trshbl = arcpy.sa.Times(rshbl, Days)

            #The second Damage mechanism due to Sheath Blight
            pshbrf = arcpy.sa.Divide(Sevl, 100)
            shbrf = arcpy.sa.Minus(1, pshbrf)
            return [rshbl, shbrf]


        #Damage mechanism due to Brown Spot
        mybsdm = 0
        bsdm = arcpy.sa.Times(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = arcpy.sa.Divide(bsdm, 100)
            pbsrf1 = arcpy.sa.Minus(1, pbsrf)
            bsrf = arcpy.sa.Power(pbsrf1, beta)
            return bsrf

        #Damage mechanism due to Bacterial Leaf Blight
        if analysisType == "Attainable Yield":
            blight = 0
            blbrf = arcpy.sa.Minus(oneRaster, blight)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLIGHT")
            for blight in BlightGDB:
                blightDate = blight[7:]
                if tempDate == blightDate:
                    pblbrf = arcpy.sa.Divide(blight, 100)
                    blbrf = arcpy.sa.Minus(1, pblbrf)
                    break
            else:
                blbrf = 1

        #Damage mechanism due to Rice Leaf Blast
        if analysisType == "Attainable Yield":
            blast = 0
            rlbrf = arcpy.sa.Minus(oneRaster, blast)

        elif analysisType == "Actual Yield":
            arcpy.AddMessage("  Calculating Reduction factor for BLAST")
            for blast in BlastGDB:
                blastDate = blast[6:]
                if tempDate == blastDate:
                    prlbrf = arcpy.sa.Divide(blast, 100)
                    prlbrf2 = arcpy.sa.Minus(1, prlbrf)
                    rlbrf = arcpy.sa.Power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1


        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        def SHR ():
            shrrf = arcpy.sa.Minus(1, shrdm)
            #shrrf = arcpy.sa.Times(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        def WH ():
            whrf = arcpy.sa.Minus(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = arcpy.sa.Minus(1, prfwd1)
            weedrf = arcpy.sa.Minus(1, rfwd)
            return weedrf

        #Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight
        def LAI ():
            arcpy.AddMessage("  Calculating Leaf Area Index")
            thisSLA = sla()
            this_shb1 = SHB()
            this_shbrf = this_shb1[1]
            this_bsrf = BS()

            if i == 1:
                pLAI = arcpy.sa.Times(thisSLA, LEAFW)
            else:
                pLAI = arcpy.sa.Times(thisSLA, thisLEAFW)

            pLAI1 = arcpy.sa.Times(blbrf, pLAI)
            pLAI2 = arcpy.sa.Times(rlbrf, pLAI1)
            pLAI3 = arcpy.sa.Times(this_shbrf, pLAI2)
            LAI = arcpy.sa.Times(this_bsrf, pLAI3)
            return LAI

        #Calculating actual radiation
        for rad in RadGDB:
            radDate = rad[3:]
            if tempDate == radDate:
                this_RAD = rad
                break


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
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt07 = (-435 * (1.7**2))+(2755 * 1.7)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where07 = "\"VALUE\" <= %d" % sumt07
            copartl = arcpy.sa.Con(co_sumt, 0.55, (arcpy.sa.Con(co_sumt, 0.45, 0, where07)), where0)
            return copartl

        #Calculating the coefficient of partitioning of roots relative to DVS
        def COPARTR ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt08 = (-435 * (1.8**2))+(2755 * 1.8)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            copartr = arcpy.sa.Con(co_sumt, 0.3, arcpy.sa.Con(co_sumt, 0.1, 0, where08), where0)
            return copartr

        #Calculating the coefficient of partitioning of panicles relative to DVS
        def COPARTP ():
            sumt075 = (-435 * (1.75**2))+(2755 * 1.75)-2320
            sumt11 = (-435 * (2.1**2))+(2755 * 2.1)-2320
            where075 = "\"VALUE\" > %d" % sumt075
            where11 = "\"VALUE\" >= %d" % sumt11
            copartp = arcpy.sa.Con(co_sumt, 1, arcpy.sa.Con(co_sumt, 0.3, 0, where075), where11)
            return copartp

        #Calculating the coefficient of partitioning of stems relative to DVS
        this_copartp = COPARTP ()
        this_copartl = COPARTL ()
        def COPARTST ():
            pr_copartst = arcpy.sa.Minus(this_copartl, this_copartp)
            copartst = arcpy.sa.Minus(1, pr_copartst)
            return copartst

        #Partitioning of assimilates
        this_copartst = COPARTST()
        this_POOL = POOL()

        def PART ():
            arcpy.AddMessage("  Calculating partitioning of assimilates")
            #Calculate the partitioning of assimilates to leaves
            this_copartr = COPARTR()
            this = arcpy.sa.Minus(1, this_copartr)
            pr_partl = arcpy.sa.Times(this_POOL, this_copartl)
            partl = arcpy.sa.Times(pr_partl, this)

            #Calculate the partitioning of assimilates to panicles
            pr_partp = arcpy.sa.Times(this_POOL, this_copartp)
            partp = arcpy.sa.Times(pr_partp, this)

            #Calculate the partitioning of assimilates to stems
            pr_partst = arcpy.sa.Times(this_POOL, this_copartst)
            partst = arcpy.sa.Times(pr_partst, this)

            #Calculate the partitioning of assimilates to roots
            partr = arcpy.sa.Times(this_POOL, this_copartr)

            return[partl, partp, partst, partr]

        #Redistribution of reserves accumulated in the stems
        def RDIST():
            arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            where1 = "\"VALUE\" < %d" % sumt1
            rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
            return rdist

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

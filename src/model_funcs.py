def calc_init_temp_sum_seed(path, IniTemp, TBASE, rt):
    """ Calculate Sum of Initial Temperature for Direct Seeded Systems

    Args:
        path (_type_): _description_
        IniTemp (_type_): _description_
        TBASE (_type_): _description_
        rt (_type_): _description_

    Returns:
        _type_: _description_
    """
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
    return rt

def calc_init_temp_sum_transplants(path, TransplantDays, IniTemp, TBASE, jt, rt, kt):
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
    return out2

def calc_sum_temp_above_TBASE(TempGDB, TBASE, out2):
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

def WH (whdm):
            whrf = arcpy.sa.Minus(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf
def WEEDS (weeddm):
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = arcpy.sa.Minus(1, prfwd1)
            weedrf = arcpy.sa.Minus(1, rfwd)
            return weedrf
        
def LAI (i, LEAFW, thisLEAFW, blbrf, rlbrf):
    """Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight"""
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

def SHR (shrdm):
    shrrf = arcpy.sa.Minus(1, shrdm)
    #shrrf = arcpy.sa.Times(pshrrf, Days)
    return shrrf

def calc_actual_radiation(RadGDB, tempDate):
        #Calculating actual radiation
        for rad in RadGDB:
            radDate = rad[3:]
            if tempDate == radDate:
                this_RAD = rad
                break
        return this_RAD

def COPARTL (co_sumt, 
             sumt00 = (-435 * (1**2))+(2755 * 1)-2320, 
             sumt07 = (-435 * (1.7**2))+(2755 * 1.7)-2320):
            #Calculating the coefficient of partitioning of leaves relative to DVS
            # sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            # sumt07 = (-435 * (1.7**2))+(2755 * 1.7)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where07 = "\"VALUE\" <= %d" % sumt07
            copartl = arcpy.sa.Con(co_sumt, 0.55, (arcpy.sa.Con(co_sumt, 0.45, 0, where07)), where0)
            return copartl

def COPARTR (co_sumt):
            #Calculating the coefficient of partitioning of roots relative to DVS
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt08 = (-435 * (1.8**2))+(2755 * 1.8)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            copartr = arcpy.sa.Con(co_sumt, 0.3, arcpy.sa.Con(co_sumt, 0.1, 0, where08), where0)
            return copartr


def COPARTP (co_sumt):
            #Calculating the coefficient of partitioning of panicles relative to DVS
            sumt075 = (-435 * (1.75**2))+(2755 * 1.75)-2320
            sumt11 = (-435 * (2.1**2))+(2755 * 2.1)-2320
            where075 = "\"VALUE\" > %d" % sumt075
            where11 = "\"VALUE\" >= %d" % sumt11
            copartp = arcpy.sa.Con(co_sumt, 1, arcpy.sa.Con(co_sumt, 0.3, 0, where075), where11)
            return copartp

        
def COPARTST(co_sumt):
    #Calculating the coefficient of partitioning of stems relative to DVS
    pr_copartst = arcpy.sa.Minus(COPARTL(co_sumt), COPARTP(co_sumt))
    copartst = arcpy.sa.Minus(1, pr_copartst)
    return copartst


def PART (co_sumt):
        #Partitioning of assimilates
        this_copartst = COPARTST(co_sumt)
        this_POOL = POOL()

        arcpy.AddMessage("  Calculating partitioning of assimilates")
        #Calculate the partitioning of assimilates to leaves
        this_copartr = COPARTR(co_sumt)
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

def RDIST(co_sumt):
#Redistribution of reserves accumulated in the stems
    arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
    sumt1 = (-435 * (2**2))+(2755 * 2)-2320
    where1 = "\"VALUE\" < %d" % sumt1
    rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
    return rdist

def sla(co_sumt):
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

def SHB (Sevl, i, LEAFW, thisLEAFW):
    #The first Damage mechanism due to Sheath Blight
    #? call rshbl, shbrf = SHB(Sevl, i, LEAFW, thisLEAFW)
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
    


def BS (beta=6.3, oneRaster):
    #! another 1* 0 raster

    #Damage mechanism due to Brown Spot
    mybsdm = 0
    bsdm = arcpy.sa.Times(oneRaster, mybsdm)

    pbsrf = arcpy.sa.Divide(bsdm, 100)
    pbsrf1 = arcpy.sa.Minus(1, pbsrf)
    bsrf = arcpy.sa.Power(pbsrf1, beta)
    return bsrf


def damage_mech_bacterial_leaf_blight(analysisType, oneRaster, BlightGDB, tempDate ):
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
        return blbrf

def damage_mech_rice_leaf_blast(analysisType, oneRaster, BlastGDB, tempDate, blast_start_date_idx=6 ):
    #Damage mechanism due to Rice Leaf Blast
    if analysisType == "Attainable Yield":
        blast = 0
        rlbrf = arcpy.sa.Minus(oneRaster, blast)

    elif analysisType == "Actual Yield":
        arcpy.AddMessage("  Calculating Reduction factor for BLAST")
        for blast in BlastGDB:
            blastDate = blast[blast_start_date_idx:]
            if tempDate == blastDate:
                prlbrf = arcpy.sa.Divide(blast, 100)
                prlbrf2 = arcpy.sa.Minus(1, prlbrf)
                rlbrf = arcpy.sa.Power(prlbrf2, 3)
                break
    else:
        rlbrf = 1
    return rlbrf

    #Calculating the Rate of Growth and the Pool of assimilates
def POOL (co_sumt,this_RAD, k=-0.6):
    arcpy.AddMessage("  Calculating Pool of assimilates")
    # sumt09 = (-435 * (1.9**2))+(2755 * 1.9)-2320
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

#Calculating the relative rate of senescence
def RRSENL (param1 = -435, param2 = 2755,param3 = 2320 ):
    arcpy.AddMessage("  Calculating the relative rate of senescence")
    #sumt107 = (-425 * (2.07**2))+(2550 * 2.07)-2125
    sumt13 = (param1 * (2.3**2))+(param2 * 2.3)-param3
    #sumt163 = (-425 * (2.63**2))+(2550 * 2.63)-2125
    where1 = "\"VALUE\" < %d" % sumt13
    #where2 = "\"VALUE\" <= %d" % sumt13
    #where3 = "\"VALUE\" <= %d" % sumt163
    rrsenl = arcpy.sa.Con(co_sumt, 0, 0.04, where1)
    return rrsenl

def RDIST(co_sumt, param1 = -300, param2 = 2000, param3 = 1700):
    #Redistribution of reserves accumulated in the stems
    arcpy.AddMessage("  Calculating the redistribution of reserves accumulated in the stems")
    sumt1 = (param1 * (2**2))+(param2 * 2)-param3
    where1 = "\"VALUE\" < %d" % sumt1
    rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
    return rdist
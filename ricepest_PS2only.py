#NAME: RICEPEST Spatial Model
#DESCRIPTION: Generic RICEPEST tool for calculating attainable and actual yields based on supplied parameters
#REQUIREMENTS: ArcGIS Spatial Analyst Extension
#DEVELOPED BY: Confidence Duku (AfricaRice), Adam Sparks (IRRI), Sander Zwart (AfricaRice
#BASED ON WORK DONE BY: Laetitia Willocquet and Serge Savary (IRRI)

#ORIGINAL VERSION of RICEPEST WITH UNMODIFIED PS3

from termios import VMIN
import rasterio
from rasterio.plot import show
from rasterio.enums import Resampling
import numpy as np
import logging
import numpy
import os
from pathlib import Path
logging.basicConfig(level=logging.INFO)



logging.info("              NAME:             RICEPEST Spatial Model")
logging.info("              DEVELOPED BY:     Confidence Duku (AfricaRice), Adam Sparks (IRRI) and Sander Zwart (AfricaRice)")
logging.info("              BASED ON WORK BY: Laetitia Willocquet and Serge Savary (IRRI)")
logging.info("              REQUIREMENTS:     ArcGIS Spatial Analyst Extension")


logging.info("\nSetting Environment Variables")

#obtaining directory for script and associated folders
# scriptPath = sys.argv[0]
# PathName = os.path.dirname(scriptPath)
PathName = Path.cwd()


#Setting environment variables;
output_dir = PathName / "Output"
yields_dir = PathName / "Yields"
# output_dir.mkdir()
# yields_dir.mkdir()

#analysisType = sys.argv[4]
analysisType = "Actual Yield"
cropEst = "Transplanted"
TransplantDays = 20
prodSituation = "PS2"

logging.info("RUNNING, " + prodSituation + " ," + analysisType)
# arcpy.CheckOutExtension("Spatial")

#Prefix the climate and disease data with the following to distinguish them
data_dir = PathName / "Data"
# RadGDB = arcpy.ListRasters("rad*", "GRID")
# TempGDB = arcpy.ListRasters("tmean*", "GRID")
# IniTemp = arcpy.ListRasters("i*", "GRID")
# BlastGDB = arcpy.ListRasters("*blast*", "TIF")
# BlightGDB = arcpy.ListRasters("*bblight*", "TIF")

###################################################
#####  New code to read in the data files  

# get raster paths
RadGDB = list(data_dir.rglob("rad*"))
RadGDB = [path for path in RadGDB if path.is_dir()]
TempGDB = list(data_dir.rglob("tmean*"))
TempGDB = [path for path in TempGDB if path.is_dir()]
IniTemp = list(data_dir.rglob("i*"))
IniTemp = [path for path in IniTemp if path.is_dir()]
IniTemp = [path for path in IniTemp if path.stem != 'info']
BlastGDB = list(data_dir.rglob("*blast*.tif",))
BlightGDB =list(data_dir.rglob("*bblight*.tif",))
all_rasters = BlightGDB + BlastGDB + IniTemp + TempGDB + RadGDB
logging.info(f" {len(all_rasters)} rasters found.")
raster_paths = all_rasters


def check_raster_crs_all_equal(raster_paths:list):
    crss = []
    for path in raster_paths:
        with rasterio.open(path) as src:
            crss.append(src.crs)
    crss = set(crss)
    logging.debug(crss)
    return crss

def find_min_raster_resolution(raster_paths:list):
    res = []
    for path in raster_paths:
        with rasterio.open(path) as src:
            res.append(src.res[0])
    return min(res)    

def get_min_intersecting_bounds(raster_paths)->[float,float,float,float]:
    """ returns the minimum bounds that intersect all rasters -> [left, bottom, right, top] """
    boundss = []
    for path in raster_paths:
        with rasterio.open(path) as src:
            boundss.append(src.bounds)
    mx = np.max(np.vstack(boundss), axis=0)
    mn = np.min(np.vstack(boundss), axis=0)
    intersecting_bounds = (mx[0], mx[1], mn[2], mn[3])
    return intersecting_bounds

def read_raster_using_bounds(raster_path:Path, 
                            bounds:[float,float,float,float], 
                            pixel_size:float,
                            resampling=Resampling.nearest)->np.ndarray:
    """
    read a rasterio window from a raster, with non-matching profile.
    read window is created from the geographic bounds and pixel size.
    """
    with rasterio.open(raster_path) as src:
        xmin, ymin, xmax, ymax = bounds
        transform = rasterio.transform.from_origin(xmin, ymax, pixel_size, pixel_size)
        raster_window = rasterio.windows.from_bounds(*bounds, transform)
        img = src.read(window=raster_window, 
                       boundless=False,         
                       resampling=resampling,
                       masked=True
                       
                       )
    return img

def write_raster_with_updated_profile(img:np.ndarray, 
                                    bounds:list, 
                                    pixel_size:float, 
                                    out_path:Path, 
                                    nodata:float=0):
    """
    save the raster with the updated profile
    profile updated to change the bounds and pixel size to match resampled data. 
    CRS is added, hard coded as WGS84
    raster saved in COG format: https://github.com/cogeotiff/rio-cogeo/blob/main/rio_cogeo/profiles.py
    """
    xmin, ymin, xmax, ymax = bounds
    transform = rasterio.transform.from_origin(xmin, ymax, pixel_size, pixel_size)
    profile = {
                "driver": "GTiff",
                "interleave": "pixel",
                "tiled": True,
                "blockxsize": 512,
                "blockysize": 512,
                "compress": "LZW",
                "width": img.shape[2],
                "height": img.shape[1],
                "count": img.shape[0],
                "dtype": img.dtype,
                "nodata": nodata,
                "transform": transform,
                "crs": "EPSG:4326"
                }
    with rasterio.open(out_path.with_suffix('.tiff'), "w", **profile) as dst:
        dst.write(img)

# ? how to handle nodata values?
# ? are nodata values set correctly? seem to be -9999.0, but metadata says {'nodata': -3.3999999521443642e+38}

raster_path = raster_paths[0]
assert check_raster_crs_all_equal(all_rasters), "all rasters must be in the same CRS, reproject before proceeding"
PIXEL_SIZE = find_min_raster_resolution(raster_paths)
BOUNDS = get_min_intersecting_bounds(raster_paths)
# img = read_raster_using_bounds(raster_path, BOUNDS, PIXEL_SIZE, resampling=Resampling.nearest)
# show(img, vmin=0)
# img.min()
# rasterio.open(raster_path).profile





# check as crs are equal
# find the intersection of the raster bounds
# find the smallest pixel size
#########################################
### Example given for PS2 only...

if prodSituation == "PS2":
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
        logging.info("Calculating Sum of Initial Temperature for Direct Seeded Systems")
        for raster in IniTemp:
            test = 0
            if rt == 0:
                out3 = raster - TBASE
                out2 = np.where(out3 < test, out3, 0)
                rt = rt + 1
                name = "out2_" + str(rt)
                write_raster_with_updated_profile(out2, BOUNDS, PIXEL_SIZE, out_path=output_dir/name)
            else:
                out1 = raster - TBASE
                out6 = np.where(out1 < test, out1, 0)
                out2 = out2 + out6
                rt = rt + 1
                name = "out2_" + str(rt)
                write_raster_with_updated_profile(out2, BOUNDS, PIXEL_SIZE, out_path=output_dir/name)
    elif cropEst == "Transplanted":
        logging.info("Calculating Sum of Initial Temperature for Transplanted Systems\n")
        for raster in IniTemp:
            test = 0
            if jt == 0:
                if kt != TransplantDays:
                    out3 = np.subtract(raster, TBASE)
                    # out4 = arcpy.sa.Con(out3, 0, out3, where1)
                    out4 = np.where(out3 < test, out3, 0)
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break
            else:
                if kt != TransplantDays:
                    out1 = np.subtract(raster, TBASE)
                    out6 = np.where(out1 < test, out1, 0)
                    out4 = out4 + out6
                    jt = jt + 1
                    kt = kt + 1
                else:
                    break

        for raster in IniTemp:
            test = 0
            if rt == 0:
                out3 = np.subtract(raster, TBASE)
                # out23 = arcpy.sa.Con(out3, 0, out3, where1)
                out23 = np.where(out3 < test, out3, 0)
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                # out23.save(output_dir / name)
                write_raster_with_updated_profile(out23, BOUNDS, PIXEL_SIZE, out_path=output_dir/name)
            else:
                out1 = np.subtract(raster, TBASE)
                out6 = np.where(out1 < test, out1, 0)
                out23 = out23 + out6
                rt = rt + 1
                ty = str(rt)
                name = "out23_" + ty
                name1 = "sumtr"
                out23.save(output_dir / name)
                write_raster_with_updated_profile(out23, BOUNDS, PIXEL_SIZE, out_path=output_dir/name)


        out22 = 0.785 * out4
        # out22.save(output_dir / "tshock")
        write_raster_with_updated_profile(out22, BOUNDS, PIXEL_SIZE, out_path=output_dir/"tshock")

        out2 = out23 - out22
        out2.save(output_dir / name1)
        write_raster_with_updated_profile(out2, BOUNDS, PIXEL_SIZE, out_path=output_dir/name1)


    #calculating sum of temp above TBASE

    i = 0
    for temp in TempGDB:
        test = 0
        where1 = "\"VALUE\" < %d" % test
        tempDate = temp[5:]
        #calculating sum of temperature
        bn = i + 15
        ft = str(bn)
        logging.info("SIMULATION AT " + ft + " DACE")
        logging.info("  Calculating Sum of Temperature above TBASE")
        if i == 0:
            out4 = np.subtract(temp, TBASE)
            out5 = arcpy.sa.Con(out4, 0, out4, where1)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(output_dir / name)
        else:
            DTemp = np.subtract(temp, TBASE)
            DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            co_sumt.save(output_dir / name)

        def sla():
            logging.info("  Calculating Specific Leaf Area")
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
            shbrf = np.subtract(1, pshbrf)
            return [rshbl, shbrf]


        #Damage mechanism due to Brown Spot
        mybsdm = 0
        bsdm = arcpy.sa.Times(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = arcpy.sa.Divide(bsdm, 100)
            pbsrf1 = np.subtract(1, pbsrf)
            bsrf = arcpy.sa.Power(pbsrf1, beta)
            return bsrf

        #Damage mechanism due to Bacterial Leaf Blight
        if analysisType == "Attainable Yield":
            blight = 0
            blbrf = np.subtract(oneRaster, blight)

        elif analysisType == "Actual Yield":
            logging.info("  Calculating Reduction factor for BLIGHT")
            for blight in BlightGDB:
                blightDate = blight[7:]
                if tempDate == blightDate:
                    pblbrf = arcpy.sa.Divide(blight, 100)
                    blbrf = np.subtract(1, pblbrf)
                    break
            else:
                blbrf = 1

        #Damage mechanism due to Rice Leaf Blast
        if analysisType == "Attainable Yield":
            blast = 0
            rlbrf = np.subtract(oneRaster, blast)

        elif analysisType == "Actual Yield":
            logging.info("  Calculating Reduction factor for BLAST")
            for blast in BlastGDB:
                blastDate = blast[6:]
                if tempDate == blastDate:
                    prlbrf = arcpy.sa.Divide(blast, 100)
                    prlbrf2 = np.subtract(1, prlbrf)
                    rlbrf = arcpy.sa.Power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = arcpy.sa.Times(oneRaster, myshrdm)

        def SHR ():
            shrrf = np.subtract(1, shrdm)
            #shrrf = arcpy.sa.Times(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = arcpy.sa.Times(oneRaster, mywhdm)

        def WH ():
            whrf = np.subtract(1, whdm)
            #whrf = arcpy.sa.Times(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = arcpy.sa.Times(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = arcpy.sa.Times(weeddm, -0.003)
            prfwd1 = arcpy.sa.Exp(prfwd)
            rfwd = np.subtract(1, prfwd1)
            weedrf = np.subtract(1, rfwd)
            return weedrf

        #Calculating LAI after treatment with Sheath Blight + Brown Spot + Bacterial leaf Blight
        def LAI ():
            logging.info("  Calculating Leaf Area Index")
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
            logging.info("  Calculating Pool of assimilates")
            where1 = "VALUE <= 0.9"
            RUE = arcpy.sa.Con(co_sumt, 1.1, 1, where1)
            this_LAI = LAI()
            pRG1 = arcpy.sa.Times(k, this_LAI)
            pRG2 = arcpy.sa.Exp(pRG1)
            pRG3 = np.subtract(1, pRG2)
            pRG4 = arcpy.sa.Times(pRG3, this_RAD)
            pRG5 = arcpy.sa.Times(pRG4, RUE)

            #Calculating the rate of growth after treatment with Weeds
            this_weedrf = WEEDS()
            RG = arcpy.sa.Times(pRG5, this_weedrf)

            #mRG = arcpy.sa.Times(RG, Days)
            ty = str(i)
            name = "pool" + ty
            RG.save(output_dir / name)
            return RG


        logging.info("  Calculating coefficients of partitioning")
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
            pr_copartst = np.subtract(this_copartl, this_copartp)
            copartst = np.subtract(1, pr_copartst)
            return copartst

        #Partitioning of assimilates
        this_copartst = COPARTST()
        this_POOL = POOL()

        def PART ():
            logging.info("  Calculating partitioning of assimilates")
            #Calculate the partitioning of assimilates to leaves
            this_copartr = COPARTR()
            this = np.subtract(1, this_copartr)
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
            logging.info("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            where1 = "\"VALUE\" < %d" % sumt1
            rdist = arcpy.sa.Con(co_sumt, 0, 4.42, where1)
            return rdist

        #Calculating the relative rate of senescence
        def RRSENL ():
            logging.info("  Calculating the relative rate of senescence")
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

        logging.info("  Calculating the dry weight of organs\n")

        #Calculate the increase in dry weight for panicles
        ppanw1 = arcpy.sa.Plus(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = arcpy.sa.Times(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = arcpy.sa.Times(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = np.subtract(this_partst, this_rdist)


        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = arcpy.sa.Times(this_rrsenl, LEAFW)
        else:
            rsenl = arcpy.sa.Times(this_rrsenl, thisLEAFW)

        pLEAFW = np.subtract(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = np.subtract(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = arcpy.sa.Plus(ppanw3, PANW)
            thisSTEMW = arcpy.sa.Plus(pstemw, STEMW)
            thisROOTW = arcpy.sa.Plus(ROOTW, this_partr)
            thisLEAFW = arcpy.sa.Plus(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(output_dir / name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            thisLEAFW.save(output_dir / name)

    output =  yields_dir + "Yield_PS2"
    thisPANW.save(output)
    logging.info("COMPLETED")



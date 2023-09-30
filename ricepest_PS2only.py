"""
Rasterio translation of the original ricepest.py script.

NOTES:
    arcpy.sa.Con is equivalent to np.where:
        where1 = "\"VALUE\" < %d" % test
        out2 = arcpy.sa.Con(out3, 0, out3, where1)

        ==>  out2 = np.where(out3 < test, out3, 0)

"""


#NAME: RICEPEST Spatial Model
#DESCRIPTION: Generic RICEPEST tool for calculating attainable and actual yields based on supplied parameters
#REQUIREMENTS: ArcGIS Spatial Analyst Extension
#DEVELOPED BY: Confidence Duku (AfricaRice), Adam Sparks (IRRI), Sander Zwart (AfricaRice
#BASED ON WORK DONE BY: Laetitia Willocquet and Serge Savary (IRRI)

#ORIGINAL VERSION of RICEPEST WITH UNMODIFIED PS3

from functools import partial
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
logging.info("              REQUIREMENTS:     tested on rasterio.__version__ '1.3.8' ")

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
TempGDB = list(data_dir.rglob("tmean*"))
RadGDB = [path for path in RadGDB if path.is_dir()]
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
    xmin, _, _, ymax = bounds
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
save_raster = partial(write_raster_with_updated_profile, bounds=BOUNDS, pixel_size=PIXEL_SIZE, out_path=output_dir)
# img = read_raster_using_bounds(raster_path, BOUNDS, PIXEL_SIZE, resampling=Resampling.nearest)
# show(img, vmin=0)
# img.min()
# rasterio.open(raster_path).profile

def lazy_load_rasters(raster_paths:list)->np.ndarray:
    for raster_path in raster_paths:
        img = read_raster_using_bounds(raster_path, 
                            BOUNDS, 
                            PIXEL_SIZE,
                            resampling=Resampling.nearest)
        yield img

TempGDB = lazy_load_rasters(TempGDB)
RadGDB = lazy_load_rasters(RadGDB)
IniTemp = lazy_load_rasters(IniTemp)
BlastGDB = lazy_load_rasters(BlastGDB)
BlightGDB = lazy_load_rasters(BlightGDB)


#########################################
### Example given for PS2 only...
# TODO: Read in the raster files...

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
# does this need as else:?
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


    logging.info("--- calculating sum of temp above TBASE ----------------")

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
            # out5 = arcpy.sa.Con(out4, 0, out4, where1)
            out5 = np.where(out4<test, out4, 0)
            co_sumt = out5 + out2
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            save_raster(img=co_sumt, out_path=output_dir / name)
        else:
            DTemp = np.subtract(temp, TBASE)
            # DTemp1 = arcpy.sa.Con(DTemp, 0, DTemp, where1)
            DTemp1 = np.where(DTemp<test, DTemp, 0)
            co_sumt = DTemp1 + co_sumt
            i = i + 1
            ty = str(i)
            name = "sumt" + ty
            save_raster(img=co_sumt, out_path=output_dir / name)

        def sla():
            #! parameterize co_sumt
            logging.info("  Calculating Specific Leaf Area")
            #calculating SLA
            #always add 1 to the value of sumt**
            #hence sumt007 with value of 0.07 becomes 1.07 in the equation
            sumt0 = (-435 * (1**2))+(2755 * 1)-2320
            sumt05 = (-435 * (1.5**2))+(2755 * 1.5)-2320
            # sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            # sumt2 = (-435 * (3**2))+(2755 * 3)-2320
            # where0 = "\"VALUE\" <= %d" % sumt0
            # where05 = "\"VALUE\" <= %d" % sumt05
            # # where1 = "\"VALUE\" <= %d" % sumt1
            # # where2 = "\"VALUE\" <= %d" % sumt2
            _SLA = np.where(co_sumt<=sumt05, 0.017, 0.03)
            SLA = np.where(co_sumt<=sumt0, _SLA, 0.037)
            # arcpy.BuildPyramids_management(SLA)
            #! I think arcpy.BuildPyramids_management is not needed, leaving without here
            return SLA

        mysevl = 0
        trysevl = co_sumt
        oneRaster = np.divide(trysevl, trysevl)
        Sevl = np.multiply(oneRaster, mysevl)

        def SHB ():
            #The first Damage mechanism due to Sheath Blight
            prshbl = np.multiply(0.00076, Sevl)
            if i == 1:
                rshbl = np.multiply(prshbl, LEAFW)
            else:
                rshbl = np.multiply(prshbl, thisLEAFW)
            #trshbl = np.multiply(rshbl, Days)

            #The second Damage mechanism due to Sheath Blight
            pshbrf = np.divide(Sevl, 100)
            shbrf = np.subtract(1, pshbrf)
            return [rshbl, shbrf]


        #Damage mechanism due to Brown Spot
        mybsdm = 0
        bsdm = np.multiply(oneRaster, mybsdm)

        def BS (beta=6.3):
            pbsrf = np.divide(bsdm, 100)
            pbsrf1 = np.subtract(1, pbsrf)
            bsrf = np.power(pbsrf1, beta)
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
                    pblbrf = np.divide(blight, 100)
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
                    prlbrf = np.divide(blast, 100)
                    prlbrf2 = np.subtract(1, prlbrf)
                    rlbrf = np.power(prlbrf2, 3)
                    break
            else:
                rlbrf = 1

        #Damage mechanism due to Sheath Rot
        myshrdm = 0
        shrdm = np.multiply(oneRaster, myshrdm)

        def SHR ():
            shrrf = np.subtract(1, shrdm)
            #shrrf = np.multiply(pshrrf, Days)
            return shrrf

        #Damage mechanism due to White Heads
        mywhdm = 0
        whdm = np.multiply(oneRaster, mywhdm)

        def WH ():
            whrf = np.subtract(1, whdm)
            #whrf = np.multiply(pwhrf, Days)
            return whrf

        #Damage mechanism due to Weeds
        myweeddm = 0
        weeddm = np.multiply(oneRaster, myweeddm)

        def WEEDS ():
            prfwd = np.multiply(weeddm, -0.003)
            prfwd1 = np.exp(prfwd)
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
                pLAI = np.multiply(thisSLA, LEAFW)
            else:
                pLAI = np.multiply(thisSLA, thisLEAFW)

            pLAI1 = np.multiply(blbrf, pLAI)
            pLAI2 = np.multiply(rlbrf, pLAI1)
            pLAI3 = np.multiply(this_shbrf, pLAI2)
            LAI = np.multiply(this_bsrf, pLAI3)
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
            # RUE = arcpy.sa.Con(co_sumt, 1.1, 1, where1)
            RUE = np.where(co_sumt <= 0.9, 1, 1.1)
            this_LAI = LAI()
            pRG1 = np.multiply(k, this_LAI)
            pRG2 = np.exp(pRG1)
            pRG3 = np.subtract(1, pRG2)
            pRG4 = np.multiply(pRG3, this_RAD)
            pRG5 = np.multiply(pRG4, RUE)

            #Calculating the rate of growth after treatment with Weeds
            this_weedrf = WEEDS()
            RG = np.multiply(pRG5, this_weedrf)

            #mRG = np.multiply(RG, Days)
            ty = str(i)
            name = "pool" + ty
            save_raster(img=RG, out_path=output_dir / name)
            return RG


        logging.info("  Calculating coefficients of partitioning")
        #Calculating the coefficient of partitioning of leaves relative to DVS
        def COPARTL ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt07 = (-435 * (1.7**2))+(2755 * 1.7)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where07 = "\"VALUE\" <= %d" % sumt07
            # _copartl = arcpy.sa.Con(co_sumt, 0.45, 0, where07)
            # copartl = arcpy.sa.Con(co_sumt, 0.55, _copartl, where0)
            _copartl = np.where(co_sumt <= sumt07, 0, 0.45)
            copartl = np.where(co_sumt <= sumt00, _copartl, 0.55)
            return copartl

        #Calculating the coefficient of partitioning of roots relative to DVS
        def COPARTR ():
            sumt00 = (-435 * (1**2))+(2755 * 1)-2320
            sumt08 = (-435 * (1.8**2))+(2755 * 1.8)-2320
            where0 = "\"VALUE\" <= %d" % sumt00
            where08 = "\"VALUE\" < %d" % sumt08
            # _copartr = arcpy.sa.Con(co_sumt, 0.1, 0, where08)
            # copartr = arcpy.sa.Con(co_sumt, 0.3, _copartr, where0)
            _copartr = np.where(co_sumt < sumt08, 0, 0.1)
            copartr = np.where(co_sumt <= sumt00, _copartr, 0.3)
            return copartr

        #Calculating the coefficient of partitioning of panicles relative to DVS
        def COPARTP ():
            sumt075 = (-435 * (1.75**2))+(2755 * 1.75)-2320
            sumt11 = (-435 * (2.1**2))+(2755 * 2.1)-2320
            where075 = "\"VALUE\" > %d" % sumt075
            where11 = "\"VALUE\" >= %d" % sumt11
            # _copartp = arcpy.sa.Con(co_sumt, 0.3, 0, where075)
            # copartp = arcpy.sa.Con(co_sumt, 1, _copartp, where11)
            _copartp = np.where(co_sumt > sumt075, 0, 0.3)
            copartp = np.where(co_sumt >= sumt11, _copartp, 1)
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
            pr_partl = np.multiply(this_POOL, this_copartl)
            partl = np.multiply(pr_partl, this)

            #Calculate the partitioning of assimilates to panicles
            pr_partp = np.multiply(this_POOL, this_copartp)
            partp = np.multiply(pr_partp, this)

            #Calculate the partitioning of assimilates to stems
            pr_partst = np.multiply(this_POOL, this_copartst)
            partst = np.multiply(pr_partst, this)

            #Calculate the partitioning of assimilates to roots
            partr = np.multiply(this_POOL, this_copartr)

            return[partl, partp, partst, partr]

        #Redistribution of reserves accumulated in the stems
        def RDIST():
            logging.info("  Calculating the redistribution of reserves accumulated in the stems")
            sumt1 = (-435 * (2**2))+(2755 * 2)-2320
            # where1 = "\"VALUE\" < %d" % sumt1
            rdist = np.where(co_sumt < sumt1, 4.42, 0)
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
            rrsenl = np.where(co_sumt <= sumt13, 0.04, 0)
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
        ppanw1 = np.add(this_rdist, this_partp)

        #Calculating the  dry weight of panicles after Sheath Rot infection
        this_shrrf = SHR()
        ppanw2 = np.multiply(ppanw1, this_shrrf)

        #Calculating the dry weight of panicles after Sheath Rot and White Head infections
        this_whrf = WH()
        ppanw3 = np.multiply(ppanw2, this_whrf)

        #Calculate the increase in dry weight for stems
        pstemw = np.subtract(this_partst, this_rdist)


        #Calculate the increase in dry weight in leaves
        if i == 1:
            rsenl = np.multiply(this_rrsenl, LEAFW)
        else:
            rsenl = np.multiply(this_rrsenl, thisLEAFW)

        pLEAFW = np.subtract(this_partl, rsenl)

        #Calculating the leaf dry weight after treatment with Sheath Blight DM 1a
        pleafw2 = np.subtract(pLEAFW, this_trshbl)

        if i == 1:
            thisPANW = np.add(ppanw3, PANW)
            thisSTEMW = np.add(pstemw, STEMW)
            thisROOTW = np.add(ROOTW, this_partr)
            thisLEAFW = np.add(LEAFW, pleafw2)
            ty = str(i)
            name = "leafw" + ty
            save_raster(img=thisLEAFW, out_path=output_dir / name)
        else:
            thisPANW = thisPANW + ppanw3
            thisSTEMW = thisSTEMW + pstemw
            thisROOTW = thisROOTW + this_partr
            thisLEAFW = thisLEAFW + pleafw2
            ty = str(i)
            name = "leafw" + ty
            save_raster(img=thisLEAFW, out_path=output_dir / name)

    output =  yields_dir + "Yield_PS2"
    save_raster(img=thisPANW, out_path=output)
    logging.info("COMPLETED")



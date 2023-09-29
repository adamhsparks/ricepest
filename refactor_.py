import arcpy
import numpy as np
import logging

def calculate_out1(raster, TBASE, where1):
    out1 = raster - TBASE
    out6 = arcpy.sa.Con(out1, 0, out1, where1)
    return out6

def calculate_out2(raster, TBASE, out2, where1):
    out1 = raster - TBASE
    out6 = arcpy.sa.Con(out1, 0, out1, where1)
    out2 = out2 + out6
    return out2

def calculate_out23(raster, TBASE, out23, where1):
    out1 = raster - TBASE
    out6 = arcpy.sa.Con(out1, 0, out1, where1)
    out23 = out23 + out6
    return out23

def calculate_out4(raster, TBASE, out4, test):
    out3 = raster - TBASE
    out4 = np.where(out3 < test, out3, 0)
    return out4

def calculate_sum_of_initial_temperature(rasters, TBASE, output_dir, prefix, test, where1):
    out = None
    for raster in rasters:
        out = calculate_out1(raster, TBASE, where1) if out is None else calculate_out2(raster, TBASE, out, where1)
        ty = str(rt)
        name = prefix + ty
        out.save(output_dir / name)
    return out

def calculate_sum_of_initial_temperature_transplanted(rasters, TBASE, output_dir, prefix, test, where1, TransplantDays):
    out1 = None
    out23 = None
    jt = 0
    kt = 0
    rt = 0
    for raster in rasters:
        if jt == 0:
            if kt != TransplantDays:
                out4 = calculate_out4(raster, TBASE, None, test)
                jt = jt + 1
                kt = kt + 1
            else:
                break
        else:
            if kt != TransplantDays:
                out4 = calculate_out4(raster, TBASE, out4, test)
                jt = jt + 1
                kt = kt + 1
            else:
                break

    for raster in rasters:
        if rt == 0:
            out23 = calculate_out23(raster, TBASE, None, where1)
            rt = rt + 1
            ty = str(rt)
            name = prefix + ty
            out23.save(output_dir / name)
        else:
            out23 = calculate_out23(raster, TBASE, out23, where1)
    return out23

# Define the input parameters
IniTemp = []
TBASE = 0
output_dir = ""
prefix = ""
test = 0
where1 = ""
cropEst = ""
TransplantDays = 0

# Call the appropriate function based on the cropEst parameter
if cropEst == "Direct Seeded":
    logging.info("Calculating Sum of Initial Temperature for Direct Seeded Systems\n")
    out = calculate_sum_of_initial_temperature(IniTemp, TBASE, output_dir, prefix, test, where1)
elif cropEst == "Transplanted":
    logging.info("Calculating Sum of Initial Temperature for Transplanted Systems\n")
    out = calculate_sum_of_initial_temperature_transplanted(IniTemp, TBASE, output_dir, prefix, test, where1, TransplantDays)
else:
    logging.error("Invalid cropEst parameter")
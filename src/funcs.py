def calc_sum_temp_above_TBASE(TempGDB, TBASE, path, out2):
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

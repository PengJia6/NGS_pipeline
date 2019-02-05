#=============================================================================
# Project : NGS_Pipleline
# Py Name: msiPipleLine
# Author : Peng Jia
# Date : 19-1-9
# Email : pengjia@stu.xjtu.edu.cn 
# Description : MSIsensor pipeline
#=============================================================================

import os
import pandas as pd
import argparse

def argumentProcress():
    """
    argument procress
    """
    global arguments, inputCase, configure
    arguments = {}
    parser = argparse.ArgumentParser(description='Generate scripts for multiple cases msisensor running! ')
    parser.add_argument('-w', '--workspace', required=True, type=str, nargs=1,
                        help="workspace of the analysis [required]")
    parser.add_argument('-i', '--caseinfo', required=True, type=str, nargs=1,
                        help="information of the case  [required]")
    parser.add_argument('-c', '--configure', required=True, type=str, nargs=1,
                        help="configure file of the pipeline [required]")
    args = parser.parse_args()
    arguments["workspace"] = args.workspace[0]
    arguments["caseinfo"] = args.caseinfo[0]
    arguments["configure"] = args.configure[0]
    ERRORSTAT = False
    if not os.path.isdir(arguments["workspace"]):
        os.system("mkdir " + arguments["workspace"])
    if len(os.listdir(arguments["workspace"])) > 0:
        print("[ERROR] The workspace should be empty! ")
        ERRORSTAT = True
    if not os.path.isfile(arguments["caseinfo"]):
        print("[ERROR] The input file is not exist in path '" + arguments["caseinfo"] + "'")
        ERRORSTAT = True
    if not os.path.isfile(arguments["configure"]):
        print("[ERROR] The configure is not exist in path '" + arguments["configure"] + "'")
        ERRORSTAT = True
    inputCase = pd.read_csv(arguments["caseinfo"], index_col=0)
    configure = pd.read_csv(arguments["configure"], index_col=0)
    if len(set(inputCase.index)) < len(inputCase.index):
        print("[ERROR] The caseinfo file '" + arguments[
            "caseinfo"] + ", may have problem, please make sure the case name unique")
    for case, info in inputCase.iterrows():
        if not os.path.isfile(info["Npath"]):
            print("[ERROR] The sequencing file '" + info["Npath"] + "' for  '" + case + "' is not exist!")
            ERRORSTAT = True
        if not os.path.isfile(info["Npath"]):
            print("[ERROR] The sequencing file '" + info["Tpath"] + "' for  '" + case + "' is not exist!")
            ERRORSTAT = True
    # if ERRORSTAT:
    #     return False
    print("[INFO] Initializing Successfully...")
    return True


def prepare():
    global arguments, inputCase
    worksapceAbs = os.path.abspath(arguments["workspace"]) + "/"
    os.system("mkdir " + worksapceAbs + "source")
    os.system("cp " + arguments["caseinfo"] + " " + worksapceAbs + "source" + "/")
    os.system("cp " + arguments["configure"] + " " + worksapceAbs + "source" + "/")
    os.system("mkdir " + worksapceAbs + "output")
    os.system("mkdir " + worksapceAbs + "script")
    arguments["output"] = worksapceAbs + "output" + "/"
    arguments["script"] = worksapceAbs + "script" + "/"
    for case in inputCase.index:
        os.system("mkdir " + arguments["output"] + case)
    return True

def generateScript():
    global configure, arguments, inputCase
    # print(configure)
    pbsDict = {}
    numOfCasePerPbs = configure.loc["numOfCasePerScript", "value"]
    casenum = 0
    for case in inputCase.index:
        pbsNum = casenum // int(numOfCasePerPbs) + 1
        if pbsNum not in pbsDict:
            pbsDict[pbsNum] = [case]
        else:
            pbsDict[pbsNum] += [case]
        casenum += 1
    for pbsid in pbsDict:
        print(pbsDict[pbsid])
        taskname = "msisensorMulti_" + str(pbsid)
        file = open(arguments["script"]  + taskname + ".sh", "w")

        file.write(
            "outputLogs="+taskname+".o\n" +
            "errorLogs="+taskname+".e\n" +
            "refMicroList="+configure.loc["refMicrosatellites","value"]+"\n"
            "thread="+configure.loc["thread","value"]+
            "echo ######################################################## >${outputLogs}\n"
            "echo case: " + " ".join(pbsDict[pbsid]) + " >>${outputLogs}\n"
            "echo ######################################################## >>${outputLogs}\n"
        )
        # print(pbsDict)
        for case in pbsDict[pbsid]:
            print(case)
            file.write(
                "echo  \n"
                "echo   >> ${outputLogs}\n"
                "echo =========================" + case + "======================== >> ${outputLogs}\n"
                "date >> ${outputLogs}\n" +
                "echo Processing " + case + "... >> ${outputLogs}\n" +
                configure.loc["msisensor","value"]+" -d ${refMicroList}"+
                   " -b ${thread}"+
                   " -n "+inputCase.loc[case,"Npath"]+
                   " -t "+inputCase.loc[case,"Tpath"]+
                   " -o "+arguments["output"]+case+"/"+case+" 2>>${errorLogs} 1>>${outputLogs} \n"
                "date >> ${outputLogs}\n" +
                "echo =========================" + case + "======================== >> ${outputLogs}\n"
                "echo   >> ${outputLogs}\n"
            )
        file.close()
        print("[INFO] Generating multiple processing script successfully!")
        return True
def main():
    if not argumentProcress():
        return -2
    if not prepare():
        return -3
    if not generateScript():
        return -4

if __name__ == "__main__":
    main()
    print()

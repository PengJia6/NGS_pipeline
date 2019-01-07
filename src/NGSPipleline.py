# =============================================================================
# Project : NGS_Pipleline
# Py Name: rawreadQC
# Author : 
# Date : 19-1-1
# Email : pengjia@stu.xjtu.edu.cn 
# Description : ''
# =============================================================================
import argparse
import os

import pandas as pd

global arguments, inputCase, configure



def argumentProcress():
    """
    argument procress
    """
    global arguments, inputCase, configure
    arguments = {}
    parser = argparse.ArgumentParser(description='DNA Next generation sequencing data process')
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

        if not os.path.isfile(info["R1"]):
            print("[ERROR] The sequencing file '" + info["R1"] + "' for  '" + case + "' is not exist!")
            ERRORSTAT = True
        if not os.path.isfile(info["R2"]):
            print("[ERROR] The sequencing file '" + info["R2"] + "' for  '" + case + "' is not exist!")
            ERRORSTAT = True

    checkList = ["fastqc", "fastp","samtools","bwa"]
    for tool in checkList:
        if tool not in configure.index:
            print("[ERROR] Please make sure the path of '"+tool+"' has been given!")
            ERRORSTAT = True
            continue
        if configure.loc[tool, "value"][-1]=="/":
            configure.loc[tool, "value"]=configure.loc[tool, "value"][:-1]
        if not os.path.isfile(configure.loc[tool,"value"]+"/"+tool):
            print("[ERROR] Please make sure '" + tool + "' is existing in path "+configure.loc[tool, "value"])
            ERRORSTAT = True
            continue
    if ERRORSTAT:
        return False
    return True


def prepare():
    global arguments, inputCase
    worksapceAbs = os.path.abspath(arguments["workspace"]) + "/"
    os.system("mkdir " + worksapceAbs + "source")
    os.system("cp " + arguments["caseinfo"] + " " + worksapceAbs + "source" + "/")
    os.system("cp " + arguments["configure"] + " " + worksapceAbs + "source" + "/")
    os.system("mkdir " + worksapceAbs + "output")
    os.system("mkdir " + worksapceAbs + "pbsfile")
    os.system("mkdir " + worksapceAbs + "pbsfile/fastqc")
    os.system("mkdir " + worksapceAbs + "pbsfile/alignment")
    os.system("mkdir " + worksapceAbs + "pbsfile/variantCalling")
    arguments["output"] = worksapceAbs + "output" + "/"
    arguments["pbs"] = worksapceAbs + "pbsfile" + "/"
    for case in inputCase.index:
        os.system("mkdir " + arguments["output"] + case)
        os.system("mkdir " + arguments["output"] + case + "/fastqc")
        os.system("mkdir " + arguments["output"] + case + "/alignment")
        os.system("mkdir " + arguments["output"] + case + "/variantCalling")
    return True




def generateFastqcPbs():
    global configure,arguments ,inputCase
    # print(configure)
    pbsDict={}
    numOfCasePerPbs=configure.loc["numOfCasePerFastqcPbs","value"]
    casenum=0
    for case in inputCase.index:
        pbsNum=casenum//int(numOfCasePerPbs)+1
        if pbsNum not in pbsDict:
            pbsDict[pbsNum]=[case]
        else:
            pbsDict[pbsNum] += [case]
        casenum += 1
    for pbsid in sorted(list(pbsDict.keys())):
        print(pbsDict[pbsid])
    for pbsid in pbsDict:
        print(pbsDict[pbsid])
        taskname = "fastqc_" + str(pbsid)
        file=open(arguments["pbs"]+"fastqc/"+taskname+".pbs","w")

        file.write(

            "#!/bin/bash\n"
            "#PBS -N "+taskname+"\n"
            "#PBS -l nodes=mu02:ppn=8\n"
            "#PBS -l walltime=9999:00:00\n"
            "#PBS -V\n"
            "cd $PBS_O_WORKDIR\n"
            
            "#######################################################\n"
                                          
            "#######################################################\n"
         
            "logs="+taskname+".logs\n"
            "REF="+str(configure.loc["ref","value"])+"\n"
            "thread="+ configure.loc["fastqcThread", "value"] + "\n"
            "echo ######################################################## >$logs\n"
            "echo case: "+" ".join(pbsDict[pbsid])+" >>$logs\n"
            "echo ######################################################## >>$logs\n"
        )
        for case in pbsDict[pbsid]:
            file.write(
                "echo   >> ${logs}\n"
                "echo ========================="+case+"======================== >> ${logs}\n"
                "date >> ${logs}\n"+
                "echo Running fastqc for "+case+" >> ${logs}\n"+
                configure.loc["fastqc","value"]+"/fastqc"
                 " -o "+arguments["output"]+case+"/fastqc/"+
                  " -t ${thread} "+
                  inputCase.loc[case,"R1"]+" "+inputCase.loc[case,"R2"]+"\n"
                "date >> ${logs}\n"+
                "echo Running fastp for "+case+" >> ${logs}\n"+
                configure.loc["fastp","value"]+"/fastp"
                 " -i " + inputCase.loc[case,"R1"]+ " -I "+ inputCase.loc[case,"R2"]+
                 " -o " + arguments["output"]+case+"/fastqc/"+case+"_HQ_R1.fastq.gz"+
                 " -O " + arguments["output"]+case+"/fastqc/"+case+"_HQ_R2.fastq.gz"+
                 " -j " + arguments["output"]+case+"/fastqc/"+case+"_fastp.json"
                 " -h " + arguments["output"]+case+"/fastqc/"+case+"_fastp.html"
                 " -w ${thread}\n" 

                                                                      
                "date >> ${logs}\n"+
                "echo ========================="+case+"======================== >> ${logs}\n"
                )
        file.close()
    return True
def generateAlignmentPbs():
    global configure, arguments, inputCase
    # print(configure)
    pbsDict = {}
    numOfCasePerPbs = configure.loc["numOfCasePerAlignmentPbs", "value"]
    casenum = 0
    for case in inputCase.index:
        pbsNum = casenum // int(numOfCasePerPbs) + 1
        if pbsNum not in pbsDict:
            pbsDict[pbsNum] = [case]
        else:
            pbsDict[pbsNum] += [case]
        casenum += 1
    for pbsid in sorted(list(pbsDict.keys())):


        taskname = "NGSAlignment_" + str(pbsid)
        file = open(arguments["pbs"] + "alignment/"+taskname + ".pbs", "w")
        # print(arguments["output"])
        file.write(

            "#!/bin/bash\n"
            "#PBS -N " + taskname + "\n"
            "#PBS -l nodes=mu02:ppn=8\n"
            "#PBS -l walltime=9999:00:00\n"
            "#PBS -V\n"
            "cd $PBS_O_WORKDIR\n"
            "#######################################################\n"
            "export BWA="+configure.loc["bwa","value"]+"/bwa\n"
            "export SAMTOOLS="+configure.loc["samtools","value"]+"/samtools\n"
            "export BAMMARKDUPLICATES="+configure.loc["bammarkduplicates","value"]+"/bammarkduplicates\n"
            "#######################################################\n"
            "logs=" + taskname + ".logs\n"
            "REF=" + str(configure.loc["ref", "value"]) + "\n"
            "thread=" +configure.loc["fastqcThread", "value"] + "\n"
            "fastqPath="+arguments["output"]+"\n"
            "echo ######################################################## >>${logs}\n"
            "echo case: " + " ".join(pbsDict[pbsid]) + " >>$logs\n"
            "echo ######################################################## >>${logs}\n"
            "caseList=("+" ".join(pbsDict[pbsid])+")\n"
            "IDList=("+" ".join(list(inputCase.loc[pbsDict[pbsid],"ID"]))+")\n"
            "SMList=("+" ".join(list(inputCase.loc[pbsDict[pbsid],"SM"]))+")\n"
            "PLList=("+" ".join(list(inputCase.loc[pbsDict[pbsid],"PL"]))+")\n"
            "LBList=("+" ".join(list(inputCase.loc[pbsDict[pbsid],"LB"]))+")\n"
            "casenum=0\n"
            "for case in ${caseList[@]}\n"
            "do\n"
            "   rg=@RG\\\\"+"tID:${IDList[$casenum]}\\\\"+"tPL:${PLList[$casenum]}\\\\"+"tSM:${SMList[$casenum]}\\\\"+"tLB:${LBList[$casenum]}"+"\n"
            "   let casenum+=1\n"
            "   echo ===========================  ${case}  ========================  >> ${logs}\n"
            "   date  >> ${logs}\n"
            "   echo Process case $(($casenum+1)) ${case} >> ${logs}\n"
            "   R1="+arguments["output"]+"${case}/fastqc/"+"${case}_HQ_R1.fastq.gz\n"
            "   R2="+arguments["output"]+"${case}/fastqc/"+"${case}_HQ_R2.fastq.gz\n"
            "   bamPATH="+arguments["output"]+"${case}/alignment/\n"
            "   echo ----------------------- bwa and sort--------------------------  >> ${logs}\n"
            "   date >> ${logs}\n"
            "   echo **${case}**  Bwa and sort ... >> ${logs}\n"
            "   ${BWA} mem -M -t ${thread} -R ${rg} ${REF} ${R1} ${R2} |${SAMTOOLS} view -@ ${thread} -Sb -|"
                                              "${SAMTOOLS} sort -@ ${thread} -o ${bamPATH}${case}_sorted.bam\n"
            "   date >> ${logs}\n"
            "   echo ----------------remove duplicated reads--------------------  >> ${logs}\n"
            "   date >> ${logs}\n"
            "   echo **${case}**  remove duplicated reads ... >> ${logs}\n"
            "   ${BAMMARKDUPLICATES} I=${bamPATH}${case}_sorted.bam O=${bamPATH}${case}_sorted_RmDup.bam "
                                              "M=${bamPATH}${case}_sorted_RG_metrics.txt "
                                              "D=${bamPATH}${case}_sorted_RG_dup.bam rmdup=0 markthreads=${thread}\n"
            "   date >> ${logs}\n"
            "   echo --------------------------bam index-------------------------  >> ${logs}\n"
            "   date >> ${logs}\n"
            "   echo **${case}**  bam index ... >> ${logs}\n"
            "   ${SAMTOOLS} index -@ ${thread} O=${bamPATH}${case}_sorted_RmDup.bam\n"
            "   date >> ${logs}\n"                                                        
            "   echo -------------- get GATK Realignment interval ---------------  >> ${logs}\n"
            "   date >> ${logs}\n"
            "   echo **${case}**  get GATK Realignment interval  ... >> ${logs}\n"
            "   java -jar "+configure.loc["javaGATK","value"]+" "
                 "-T RealignerTargetCreator " 
                 "-R ${REF} "
                 "-I ${bamPATH}${case}_sorted_RmDup.bam " 
                 "-known "+configure.loc["1KGP3Indels","value"]+" " 
                 "-known "+configure.loc["1KGGoldIndels","value"]+" "
                 "-o ${bamPATH}${case}_IndelRealigner.intervals "
                 "-nct ${thread}\n"   
            "   date >> ${logs}\n"                                                        
            "   echo -------------- GATK Realignment  ---------------  >> ${logs}\n"
            "   date >> ${logs}\n"
            "   echo **${case}** GATK Realignment  ... >> ${logs}\n"
            "   java -jar "+configure.loc["javaGATK","value"]+" "
                 "-T IndelRealigner " 
                 "-R ${REF} "
                 "-I ${bamPATH}${case}_sorted_RmDup.bam " 
                 "-known "+configure.loc["1KGP3Indels","value"]+" " 
                 "-known "+configure.loc["1KGGoldIndels","value"]+" "
                 "-o ${bamPATH}${case}_sorted_RmDup_realign.bam "
                 "--targetIntervals ${bamPATH}${case}_IndelRealigner.intervals "
                 "-nct ${thread}\n"
            "   date >> ${logs}\n"  
            "   echo ------------------ get BRSR table  ---------------------  >> ${logs}\n"
            "   date >> ${logs}\n"
            "   echo **${case}** get BRSR table  ... >> ${logs}\n"
            "   java -jar "+configure.loc["javaGATK","value"]+" "
                 "-T BaseRecalibrator " 
                 "-R ${REF} "
                 "-I ${bamPATH}${case}_sorted_RmDup_realign.bam " 
                 "-knownSites "+configure.loc["1KGP3Indels","value"]+" " 
                 "-knownSites "+configure.loc["1KGGoldIndels","value"]+" "
                 "-knownSites "+configure.loc["dnsnp","value"]+" "
                 "-o ${bamPATH}${case}_BQSR.table "
                 "-nct ${thread}\n"
            "   date >> ${logs}\n"
            "   echo ------------------- BRSR  ------------------------------  >> ${logs}\n"
            "   date >> ${logs}\n"
            "   echo **${case}** get BRSR table  ... >> ${logs}\n"
            "   java -jar "+configure.loc["javaGATK","value"]+" "
                 "-T PrintReads " 
                 "-R ${REF} "
                 "-I ${bamPATH}${case}_sorted_RmDup_realign.bam " 
                 "-BQSR ${bamPATH}${case}_BQSR.table "
                 "-o ${bamPATH}${case}_sorted_RmDup_realign_BQSR.bam "
                 "-nct ${thread}\n"                
            "   date >> ${logs}\n"            
            "   echo --------------------------------------------------------  >> ${logs}\n"                     
            "   date >> ${logs}\n"
            "   echo =============================================================== >> ${logs}\n"                                   
            "done\n"
                                                  
                                                  
                                                  
#                                                   
#                                                   " date >> ${logs}"
#     bwa mem -M -t ${thread} ${REF} ${fastqPath}${sample}_R1.fastq.gz ${fastqPath}${sample}_R2.fastq.gz |samtools view -@ ${thread} -Sb -|samtools sort -@ ${thread} -o ${sample}_sorted.bam
#     date >> ${logs}
#     echo ----------------------------------------------  >> ${logs}
#     echo -------------------add read group-------------  >> ${logs}
#     date >> ${logs}
#     samtools addreplacerg -@ ${thread} -r 'ID:XJTU' -r 'SM:RPG' -r 'LB:DNA' -r 'PL:ILM' ${sample}_sorted.bam |samtools view -@ 60 -Sb -o ${sample}_sorted_RG.bam -
#     date >> ${logs}
#     echo ----------------------------------------------  >> ${logs}
#     echo -----------------remove duplicates------------  >> ${logs}
#     date >> ${logs}
#     bammarkduplicates I=${sample}_sorted_RG.bam O=${sample}_sorted_RG_RmDup.bam M=${sample}_sorted_RG_metrics.txt D=${sample}_sorted_RG_dup.bam rmdup=1 markthreads=${thread}
#     date >> ${logs}
# "
        )

    return True
def main():
    if not argumentProcress():
        return -2
    if not prepare():
        return -3
    if not generateFastqcPbs():
        return -4
    if not generateAlignmentPbs():
        return -5
if __name__ == "__main__":
    main()
    print()
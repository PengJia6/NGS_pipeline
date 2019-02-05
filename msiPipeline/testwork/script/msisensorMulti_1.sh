outputLogs=msisensorMulti_1.o
errorLogs=msisensorMulti_1.e
refMicroList=/home/REF/GRCh38.d1.vd1.fa.microsatellites
thread=6echo ######################################################## >${outputLogs}
echo case: a b c >>${outputLogs}
echo ######################################################## >>${outputLogs}
echo  
echo   >> ${outputLogs}
echo =========================a======================== >> ${outputLogs}
date >> ${outputLogs}
echo Processing a... >> ${outputLogs}
/home/pengjia/anaconda3/bin/msisensor -d ${refMicroList} -b ${thread} -n a_R1 -t aR2 -o /mnt/project/MyPipline/github_NGS_Pipeline/NGS_pipeline/msiPipeline/testwork/output/a/a 2>>${errorLogs} 1>>${outputLogs} 
date >> ${outputLogs}
echo =========================a======================== >> ${outputLogs}
echo   >> ${outputLogs}
echo  
echo   >> ${outputLogs}
echo =========================b======================== >> ${outputLogs}
date >> ${outputLogs}
echo Processing b... >> ${outputLogs}
/home/pengjia/anaconda3/bin/msisensor -d ${refMicroList} -b ${thread} -n b_R1 -t bR2 -o /mnt/project/MyPipline/github_NGS_Pipeline/NGS_pipeline/msiPipeline/testwork/output/b/b 2>>${errorLogs} 1>>${outputLogs} 
date >> ${outputLogs}
echo =========================b======================== >> ${outputLogs}
echo   >> ${outputLogs}
echo  
echo   >> ${outputLogs}
echo =========================c======================== >> ${outputLogs}
date >> ${outputLogs}
echo Processing c... >> ${outputLogs}
/home/pengjia/anaconda3/bin/msisensor -d ${refMicroList} -b ${thread} -n cN -t cT -o /mnt/project/MyPipline/github_NGS_Pipeline/NGS_pipeline/msiPipeline/testwork/output/c/c 2>>${errorLogs} 1>>${outputLogs} 
date >> ${outputLogs}
echo =========================c======================== >> ${outputLogs}
echo   >> ${outputLogs}

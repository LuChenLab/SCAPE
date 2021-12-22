# conda create -n polyApipe_env  --override-channels -c bioconda -c conda-forge -c anaconda umi_tools=1.0.0-0 pysam samtools subread 
conda activate polyApipe_env

/home/zhouran/data/proj/2021-0918-apa_evaluation/04.evaluation_polyA/polyApipe/code/polyApipe-master/polyApipe.py -i \
/home/zhouran/data/proj/2021-0918-apa_evaluation/03.simulation_data_with_polyA/align/simu.Aligned.sortedByCoord.out.addtag.bam -o simu

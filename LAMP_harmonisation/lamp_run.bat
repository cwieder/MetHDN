@echo off

set id=%MTBLS3852
set ion_mode=%pos

lamp cli ^
  --sep "tab" ^
  --input-data "D:\MARIANA2\HDN\MetHDN\LAMP_harmonisation\%id%\input\xcms_data.tsv" ^
  --col-idx "1, 3, 6, 12" ^
  --add-path "" ^
  --ref-path "D:\MARIANA2\HDN\MetHDN\LAMP_harmonisation\lam_ebi\lamp_ebi\test-data\ref_all_v7_pos.tsv" ^
  --ion-mode "%ion_mode%" ^
  --thres-rt "1.0" ^
  --thres-corr "0.5" ^
  --thres-pval "0.05" ^
  --method "pearson" ^
  --positive ^
  --ppm "5.0" ^
  --save-mr ^
  --sr-out "D:\MARIANA2\HDN\MetHDN\LAMP_harmonisation\%id%\res\%id%_LAMP.tsv" ^
  --mr-out "D:\MARIANA2\HDN\MetHDN\LAMP_harmonisation\%id%\res\%id%_LAMP_mr.tsv" ^
  --db-out "D:\MARIANA2\HDN\MetHDN\LAMP_harmonisation\%id%\res\%id%_LAMP_db.sqlite" ^

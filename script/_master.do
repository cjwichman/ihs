
*****************************************************************
* Replication file for:
* "Elasticities and the inverse hyperbolic sine transformation"
* by Marc F. Bellemare & Casey J. Wichman (2019)
*****************************************************************



clear all
drop _all
scalar drop _all



********  change this directory  *********
global dir = "/Users/wichman/Dropbox/IHS" 
******************************************



****   Main analysis   ****
cd $dir
run "script/makeTable1.do"
run "script/makeTable2.do"
run "script/makeTable3.do"
***************************

* supplementary code for Appendix B
doedit "script/IHSElasticities.do"

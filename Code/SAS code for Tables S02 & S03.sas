options ls=200 ps=60 nodate nonumber msglevel=I validvarname=any 
	validmemname=extend;
libname xl "/home/liux3/Manuscript/26/Revised Manuscript";
****** Read excel sheets into SAS;
libname myexcel xlsx 
	"/home/liux3/Manuscript/26/Revised Manuscript/Data/Final new data.xlsx";

/*proc contents data=myexcel._all_;
run;*/
****** Create macro to calculate adjusted aplph values and to conduct MANOVA analysis;
*** Basic macro function;

%macro Liu(fig, vars);
	data fig_1;
		set myexcel.&fig;

		if diet=' ' then
			delete;
	run;

	proc princomp data=fig_1 out=fig_1pc outstat=stats noprint;
		var &vars;
	run;

	*** Get eigenvalues from data named stats from PC analysis;

	data eigenvalues;
		set stats;
		where _type_='EIGENVAL';
		drop _NAME_;
		* drop this temp variable because of missing values to avoid listed issues in proc transpose below;
	run;

	proc transpose data=eigenvalues out=eigen_long(rename=(COL1=eigenvalue));
		var &vars;
	run;

	*** Calculate adjusted alpha using Cheverud-Repsilber eigenvalue correction;

	data alpha_correction;
		length fig_id $25;
		set eigen_long end=last;
		retain Fig_id M top_pc M_eff Original_Alpha Adjusted_Alpha1 Adjusted_Alpha2;
		retain sum_excess 0 M 0 top_pc 0;
		*** Count total variables;
		M + 1;
		*** Sum up the variance excess for all eigenvalues greater than 1;

		if eigenvalue > 1 then
			do;
				sum_excess=sum_excess + (eigenvalue - 1);
				top_pc+1;
			end;

		if last then
			do;
				M_eff=M - sum_excess;
				Original_Alpha=0.05;
				*** Apply the Cheverud-Repsilber eigenvalue correction;
				Adjusted_Alpha1=Original_Alpha / top_pc;
				Adjusted_Alpha2=Original_Alpha / M_eff;				
				*** Add a new variable here;
				Fig_id="%sysfunc(dequote(&fig))";
				output;
			end;
		keep Fig_id M M_eff Original_Alpha Adjusted_Alpha1 Adjusted_Alpha2 top_pc;
		*proc print;
	run;

	****** Conduct multivariate ANOVA (MANOVA) analysis;
	ods trace on;
	ods select none;
	ods output multstat=manova_p1;
	%let clean_fig = %sysfunc(substr(%sysfunc(dequote(&fig)), 1, 6));

	%if &clean_fig=%str(Fig 6 ) %then
		%do;

			proc glm data=fig_1;
				class diet(ref="Control");
				model &vars=diet;
				manova h=_all_/printe printh;
				* Request multivariate tests for all effects;
				run;
			quit;

			data manova_p(rename=(hypothesis=factor));
				retain Fig_id;
				set manova_p1;

				if statistic^="Wilks' Lambda" then
					delete;
				pvalue=ProbF;
				fig_id="%sysfunc(dequote(&fig))";
				format pvalue best10.;
				keep Fig_id hypothesis pvalue;
			run;

		%end;
	%else %if &clean_fig=Fig S1 %then
		%do;

			proc glm data=fig_1;
				class diet(ref="Palm");
				* Your categorical independent variable;
				model &vars=diet;
				manova h=_all_/printe printh;
				* Request multivariate tests for all effects;
				run;
			quit;

			data manova_p(rename=(hypothesis=factor));
				retain Fig_id;
				set manova_p1;

				if statistic^="Wilks' Lambda" then
					delete;
				pvalue=ProbF;
				fig_id="%sysfunc(dequote(&fig))";
				format pvalue best10.;
				keep Fig_id hypothesis pvalue;
			run;

		%end;
	%else
		%do;

			proc glm data=fig_1;
				class 'germ status'n(ref='MPF') diet(ref='Palm');
				* Your categorical independent variable;
				model &vars='germ status'n diet 'germ status'n * diet;
				manova h=_all_/printe printh;
				* Request multivariate tests for all effects;
				run;
			quit;

			ods output multstat=manova_p2;

			proc glm data=fig_1;
				class 'One way group'n(ref='MPF Palm');
				* Your categorical independent variable;
				model &vars='One way group'n;
				manova h=_all_/printe printh;
				* Request multivariate tests for all effects;
				run;
			quit;

			data manova_p(rename=(hypothesis=factor));
				retain Fig_id;
				set manova_p1 manova_p2;

				if hypothesis="ONE WAY GROUP" then
					hypothesis="OVERALL";

				if statistic^="Wilks' Lambda" then
					delete;
				pvalue=ProbF;
				fig_id="%sysfunc(dequote(&fig))";
				format pvalue best10.;
				keep Fig_id hypothesis pvalue;
				*proc print;
			run;

		%end;
	ods select all;
	ods trace off;

	proc transpose data=manova_p out=wide_manova_p;
		by fig_id;
		id factor;
		var pvalue;
	run;

	data exact_p;
		length fig_id $25;
		set wide_manova_p;
		drop _name_ _label_;
		*proc print;
	run;

%mend Liu;

*** Analysis of one class of variables in each figure;

%macro one_class(onefig);
%global fig;
	%local num_fid i this_id;
	%let num_fid=%sysfunc(countw(&onefig, %str(|)));

	%do i=1 %to &num_fid;
		%let this_id=%scan(&onefig, &i, %str(|));

		%if &i=1 %then %do;
			%let fig=&this_id;
			           %end;
		%else
			%do;
				%Liu(&fig, &this_id);
			%end;
	%end;
%mend one_class;
/*
%one_class('Fig 1D'n|'Liver Weight (g)'n);

%one_class('FIGURE 1 E-J'n|'Plasma Active GLP-1'n 'Plasma Glucagon'n 
	'Plasma PYY'n 'Plasma Insulin'n 'Plasma MCP-1'n 'Plasma Leptin'n);
	
%one_class('Fig 1 K-L revised'n|'gWAT Weight'n 'Fast Glucose'n);

%one_class('Fig 1 M revised'n|'Steatosis %'n);

%one_class('Fig3 A-H revised'n|'Plasma FA 22:5'n 'Plasma TG 18:1_20:4_22:6'n 'Plasma TG 16:0_18:1_18:1 isomer'n 'Plasma KDdiA-PC'n 
	'Plasma Cer d 18:1_24:1'n 'Plasma Cer.d.18:1-22:0'n 'Plasma Acylcarnitine 16:1'n 'Plasma Cholesterol'n);
	
%one_class('Fig3 I-P'n|'Liver PE 38:4'n 'Liver PE 38:5'n 'Liver PC38:5 IsomerA'n 'Liver PC40:6'n 'Liver PI 38:3'n 'Liver Cer d42-2'n 
	'Liver LPE-N  36:4'n 'Liver TG O57:3'n);
	
%one_class('Fig 4 A-H'n|'LIVER aHYA'n 'LIVER aKETOD'n 'LIVER CLA1'n 'LIVER HYE'n 'LIVER gKETOB'n 'LIVER gHYA'n 'LIVER HYARA'n 'LIVER 12-OH'n);

%one_class('Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n 'LIVER Protectin dx'n 'LIVER LIPOXIN'n 'LIVER RESOLVIN D5'n 'LIVER 20-HETE'n 'LIVER N-DGLA EA'n 
	'LIVER 14,15 DIHETE'n 'LIVER 12,13 DIHOME'n);
	
%one_class('Fig 6 B-L revised'n|'Liver Ricinoleic Acid'n 'Liver 12-Hydroxystearic Acid'n 'Liver 12-Oxostearic Acid'n 'Body Weight (g)'n 
	'Liver Weight (g)'n 'gWAT Weight (g)'n 'Necropsy Blood Glucose (mg/dL)'n 'Liver TG 14:0_16:0_16:0'n 'Liver PC 20:1_14:1'n 
	'Liver 20-Hydroxydocosa-hexaenoic'n 'Liver Cer d38:1'n);
	
%one_class('Fig 6 P-S revised'n|'Rsph9'n 'Ciart'n 'Slc01a1'n 'Cyp4a14'n);

%one_class('Fig S1 revised'n|'Diet TG 12:0_12:0_14:0'n 'Diet  TG 18:1_18:1_18:1'n 'Diet TG 18:1_18:2_18:2'n 'Diet TG 18:2_20:3_22:0'n 
	'Diet  TG 18:3_18:3_18:3'n 'Diet TG 12:0_18:2_22:6'n 'Diet TG 16:0_20:4_20:5 minor'n);
	
%one_class('Fig S2 revised'n|'Bacteroidetes'n 'Deferribacteres'n 'Verrucomicrobia'n 'Desulfobacterota'n 'Delta and gammaproteobacteria'n 
	'Firmicutes'n);
	
%one_class('Fig S3 revised done'n|'Hepatic Srebf1'n 'Hepatic ACC1'n 'Hepatic FASN'n 'Hepatic CD68'n);
%one_class('Fig S3 revised done'n|'gWAT SREBF1'n 'gWAT ACC1'n 'gWAT FASN'n 'gWAT CD68'n);
%one_class('Fig S3 revised done'n|'COLON SREBF1'n 'COLON ACC1'n 'COLON FASN'n 'COLON CD68'n);
	
%one_class('Fig S5 A-H revised done'n|'PLASMA HYA'n 'PLASMA HYB'n 'PLASMA aHYA+aHYC'n 'PLASMA CALA2+CGLA1'n 'PLASMA CGLA2'n 'PLASMA gHYD'n 
	'PLASMA 12-HAS'n 'PLASMA 12-OXO'n);
	
%one_class('Fig S5 I-P revised done'n|'ILEUM HYA'n 'ILEUM KETOA'n 'ILEUM aHYA+aHYC'n 'ILEUM aKETOD'n 'ILEUM gHYD'n 'ILEUM gLA+CALA3'n 
	'ILEUM 12-HAS'n 'ILEUM 12-OXO'n);
	
%one_class('FigS6 revised done'n|'Leukotriene B4'n '11,12-DiHETre'n '9.10-DiHOME'n '14,15-DiHETre'n '16,17-DeHDoPE'n '12-HODE'n '8,9-DiHETre'n 
	'11,12-DiHETE'n '5-HEPE'n '10-HODE'n '15-deoxy-PGJ2'n '9-HOTE'n '8,9-DiHETE'n '17,18-DiHETE'n '8,15-DiHETE'n 'Prostaglandin E2'n 'N-C22:6n3_EA'n '15-HETE'n 
	'13-KODE'n '9-HODE'n);
	
%one_class('Fig S9 revised done 'n|'Glycine C18:1'n 'Glycine C12:2'n 
	'Methionine C2:0'n 'Lysine C18:1'n 'Tyrosine C 18:2'n 
	'Alanine/Sarcosine - C18:2'n);
*/
	

*****************************************************************************************************	
********************** Analysis of all classes of variables for all figures *************************
*****************************************************************************************************;

%macro all_class(allfig);
%local num_fad a this_ad;

%let num_fad=%sysfunc(countw(&allfig, %str(~)));
%do a=1 %to &num_fad;
    %let this_ad=%scan(&allfig, &a, %str(~));
    %one_class(&this_ad);
    %if &a=1 %then %do;
         data alpha_correction_all;
         set alpha_correction;
         run;
             
         data global_p_all;
         set exact_p;
         run;    
                 %end;
    %else %do;
         data alpha_correction_all;
         set alpha_correction_all alpha_correction;
         run;
             
         data global_p_all;
         set global_p_all exact_p;          
         run;   
          %end;
%end;

%mend all_class;

	
%all_class('Fig 1 D'n|'Liver Weight (g)'n~

'Fig 1 E-J'n|'Plasma Active GLP-1'n 'Plasma Glucagon'n 'Plasma PYY'n 'Plasma Insulin'n 'Plasma MCP-1'n 'Plasma Leptin'n~

'Fig 1 K-L revised'n|'gWAT Weight'n 'Fast Glucose'n~

'Fig 1 M revised'n|'Steatosis %'n~

'Fig 3 A-H revised'n|'Plasma FA 22:5'n 'Plasma TG 18:1_20:4_22:6'n 'Plasma TG 16:0_18:1_18:1 isomer'n 'Plasma KDdiA-PC'n 
	'Plasma Cer d 18:1_24:1'n 'Plasma Cer.d.18:1-22:0'n 'Plasma Acylcarnitine 16:1'n 'Plasma Cholesterol'n~

'Fig 3 I-P'n|'Liver PE 38:4'n 'Liver PE 38:5'n 'Liver PC38:5 IsomerA'n 
	'Liver PC40:6'n 'Liver PI 38:3'n 'Liver Cer d42-2'n 'Liver LPE-N  36:4'n 'Liver TG O57:3'n~

'Fig 4 A-H'n|'LIVER aHYA'n 'LIVER aKETOD'n 'LIVER CLA1'n 'LIVER HYE'n 'LIVER gKETOB'n 'LIVER gHYA'n 'LIVER HYARA'n 'LIVER 12-OH'n~

'Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n 'LIVER Protectin dx'n 'LIVER LIPOXIN'n 
	'LIVER RESOLVIN D5'n 'LIVER 20-HETE'n 'LIVER N-DGLA EA'n 'LIVER 14,15 DIHETE'n 'LIVER 12,13 DIHOME'n~
	
'Fig S2 revised'n|'Bacteroidetes'n 'Deferribacteres'n 'Verrucomicrobia'n 'Desulfobacterota'n 'Delta and gammaproteobacteria'n 'Firmicutes'n~

'Fig S3 revised done'n|'Hepatic Srebf1'n 'Hepatic ACC1'n 'Hepatic FASN'n 'Hepatic CD68'n~
'Fig S3 revised done'n|'gWAT SREBF1'n 'gWAT ACC1'n 'gWAT FASN'n 'gWAT CD68'n~
'Fig S3 revised done'n|'COLON SREBF1'n 'COLON ACC1'n 'COLON FASN'n 'COLON CD68'n~
 
'Fig S5 A-H revised done'n|'PLASMA HYA'n 'PLASMA HYB'n 'PLASMA aHYA+aHYC'n 'PLASMA CALA2+CGLA1'n 'PLASMA CGLA2'n 
         'PLASMA gHYD'n 'PLASMA 12-HAS'n 'PLASMA 12-OXO'n~
         
'Fig S5 I-P revised done'n|'ILEUM HYA'n 'ILEUM KETOA'n 'ILEUM aHYA+aHYC'n 
	'ILEUM aKETOD'n 'ILEUM gHYD'n 'ILEUM gLA+CALA3'n 'ILEUM 12-HAS'n 'ILEUM 12-OXO'n~ 

'Fig S6 revised done'n|'Leukotriene B4'n '11,12-DiHETre'n '9.10-DiHOME'n '14,15-DiHETre'n '16,17-DeHDoPE'n '12-HODE'n '8,9-DiHETre'n '11,12-DiHETE'n 
	'5-HEPE'n '10-HODE'n '15-deoxy-PGJ2'n '9-HOTE'n '8,9-DiHETE'n '17,18-DiHETE'n '8,15-DiHETE'n 'Prostaglandin E2'n 'N-C22:6n3_EA'n '15-HETE'n 
	'13-KODE'n '9-HODE'n~

'Fig S9 revised done 'n|'Glycine C18:1'n 'Glycine C12:2'n 'Methionine C2:0'n 
	'Lysine C18:1'n 'Tyrosine C 18:2'n 'Alanine/Sarcosine - C18:2'n);
	
data alpha_correction_all1;
    set alpha_correction_all;
run;

data global_p_all1;
    set global_p_all;
run;

%all_class('Fig 6 B-L revised'n|'Liver Ricinoleic Acid'n 'Liver 12-Hydroxystearic Acid'n 'Liver 12-Oxostearic Acid'n 'Body Weight (g)'n 
	'Liver Weight (g)'n 'gWAT Weight (g)'n 'Necropsy Blood Glucose (mg/dL)'n 'Liver TG 14:0_16:0_16:0'n 'Liver PC 20:1_14:1'n 
	'Liver 20-Hydroxydocosa-hexaenoic'n 'Liver Cer d38:1'n~
	
'Fig 6 P-S revised'n|'Rsph9'n 'Ciart'n 'Slc01a1'n 'Cyp4a14'n~

'Fig S1 revised'n|'Diet TG 12:0_12:0_14:0'n 'Diet  TG 18:1_18:1_18:1'n 'Diet TG 18:1_18:2_18:2'n 'Diet TG 18:2_20:3_22:0'n 
	'Diet  TG 18:3_18:3_18:3'n 'Diet TG 12:0_18:2_22:6'n 'Diet TG 16:0_20:4_20:5 minor'n);
	
data alpha_correction_all2;
    set alpha_correction_all;
run;
data global_p_all2;
    set global_p_all;
    OVERALL=diet;
run;



data alpha_correction_all;
    set alpha_correction_all1 alpha_correction_all2;
    if top_pc=0 then 
         do; 
             top_pc=1; Adjusted_Alpha1=Original_Alpha;
         end;
    proc sort; by fig_id;
run;
data global_p_all;
    set global_p_all1 global_p_all2;
    proc sort; by fig_id;
run;	

*ods rtf close;

proc export data=alpha_correction_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/alpha_correction_all.csv" 
		dbms=csv replace;
run;

proc export data=global_p_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/global_p_all.csv" 
		dbms=csv replace;
run;





****** Different way to run analysis of one class of variables in each figure;

/*%macro liu(fig, vars);
data fig_1; set myexcel.&fig;
if diet=' ' then delete;
run;

proc princomp data=fig_1 out=fig_1pc outstat=stats noprint;
var &vars;
run;

 *** Get eigenvalues from data named stats from PC analysis;
data eigenvalues;
set stats;
where _type_='EIGENVAL';
drop _NAME_; * drop this temp variable because of missing values to avoid listed issues in proc transpose below;
run;

proc transpose data=eigenvalues out=eigen_long(rename=(COL1=eigenvalue));
var &vars;
run;


 *** Calculate adjusted alpha using Cheverud-Repsilber eigenvalue correction;
data alpha_correction;
set eigen_long end=last;
retain Fig_id M top_pc M_eff Original_Alpha Adjusted_Alpha;
retain sum_excess 0 M 0 top_pc 0;

 *** Count total variables;
M + 1;

 *** Sum up the variance excess for all eigenvalues greater than 1;
if eigenvalue > 1 then do;
sum_excess = sum_excess + (eigenvalue - 1);
top_pc+1;
end;

if last then do;
M_eff = M - sum_excess;
Original_Alpha = 0.05;

 *** Apply the Cheverud-Repsilber eigenvalue correction;
Adjusted_Alpha = Original_Alpha / M_eff;

 *** Add a new variable here;
Fig_id="%sysfunc(dequote(&fig))";

output;
end;
keep Fig_id M M_eff Original_Alpha Adjusted_Alpha top_pc;
proc print;
run;


 ****** Conduct multivariate ANOVA (MANOVA) analysis;
ods trace on;
ods select none;
ods output multstat=manova_p1;

proc glm data=fig_1;
class 'germ status'n(ref='MPF') diet(ref='Palm'); * Your categorical independent variable;
model &vars = 'germ status'n diet 'germ status'n * diet;
manova h=_all_/printe printh; * Request multivariate tests for all effects;
run;
quit;

ods output multstat=manova_p2;
proc glm data=fig_1;
class 'One way group'n(ref='MPF Palm'); * Your categorical independent variable;
model &vars = 'One way group'n;
manova h=_all_/printe printh; * Request multivariate tests for all effects;
run;
quit;

ods select all;
ods trace off;

data manova_p(rename=(hypothesis=factor));
retain Fig_id;
set manova_p1 manova_p2;
if hypothesis="ONE WAY GROUP" then hypothesis="OVERALL";
if statistic^="Wilks' Lambda" then delete;
pvalue=ProbF;
fig_id="%sysfunc(dequote(&fig))";
format pvalue best10.;
keep Fig_id hypothesis pvalue;
 *proc print;
run;

proc transpose data=manova_p out=wide_manova_p;
by fig_id;
id factor;
var pvalue;
run;

data exact_p;
set wide_manova_p;
drop _name_ _label_;
proc print;
run;

%mend liu;

%liu('FIGURE 1 E-J'n, 'Plasma Active GLP-1'n 'Plasma Glucagon'n 'Plasma PYY'n 'Plasma Insulin'n 'Plasma MCP-1'n 'Plasma Leptin'n);


%liu('Fig S1 A-L'n, 'Hepatic Srebf1'n 'Hepatic ACC1'n 'Hepatic FASN'n 'Hepatic CD68'n);
%liu('Fig S1 A-L'n, 'gWAT SREBF1'n 'gWAT ACC1'n 'gWAT FASN'n 'gWAT CD68'n);
%liu('Fig S1 A-L'n, 'COLON SREBF1'n 'COLON ACC1'n 'COLON FASN'n 'COLON CD68'n);


%liu('Fig 4 A-H'n, 'LIVER aHYA'n 'LIVER aKETOD'n 'LIVER CLA1'n 'LIVER HYE'n	'LIVER gKETOB'n	'LIVER gHYA'n 'LIVER HYARA'n 'LIVER 12-OH'n);
%liu('Fig 4 I-P'n, 'LIVER 18-RESOLVIN E2'n 'LIVER Protectin dx'n 'LIVER LIPOXIN'n 'LIVER RESOLVIN D5'n 'LIVER 20-HETE'n	'LIVER N-DGLA EA'n 'LIVER 14,15 DIHETE'n
'LIVER 12,13 DIHOME'n);


%liu('Fig S4 A-H'n, 'PLASMA HYA'n 'PLASMA HYB'n	'PLASMA aHYA+aHYC'n	'PLASMA CALA2+CGLA1'n 'PLASMA CGLA2'n 'PLASMA gHYD'n 'PLASMA 12-OH'n 'PLASMA 12-OXO'n);
%liu('Fig S4 I-P'n, 'ILEUM HYA'n 'ILEUM KETOA'n	'ILEUM aHYA+aHYC'n 'ILEUM aKETOD'n 'ILEUM gHYD'n 'ILEUM gLA+CALA3'n	'ILEUM 12-OH'n 'ILEUM 12-OXO'n);


%liu('FigS5 liver oxylipins'n, 'Leukotriene B4'n '11,12-DiHETre'n '9.10-DiHOME'n '14,15-DiHETre'n '16,17-DeHDoPE'n '12-HODE'n '8,9-DiHETre'n
'11,12-DiHETE'n	'5-HEPE'n '10-HODE'n '15-deoxy-PGJ2'n '9-HOTE'n	'8,9-DiHETE'n '17,18-DiHETE'n '8,15-DiHETE'n 'Prostaglandin E2'n 'N-C22:6n3_EA'n
'15-HETE'n	'13-KODE'n	'9-HODE'n);


%liu('Fig3 plasma FA'n, '18:0 (SA)'n '18:2 (LA)'n '18:3 Gamma Linolenic Acid'n '18:3 Alpha Linolenic Acid'n	'18:4 Stearidonic Acid'n
'20:4 Arachidonic Acid'n '20:5 EPA'n '22:6 Docosahexaenoic Acid'n '16:0 PA'n '16:1 Palmitoleic Acid'n '20:3 Trienoic Acid'n '22:5 Docosapentaenoic'n);


%liu('Fig3 liver lipidomics'n, 'PE 38:4'n   'PE 38:5'n 'PC38:5 IsomerA'n 'PC40:6'n 'PI 38:3'n 'Cer d42-2'n 'LPE-N'n	'TG O57:3'n	'FA 20:3'n
'FA 22:4'n 'DG16:0_18:1'n 'DG16:0_18:2'n 'FA 18:3'n	'PC38:4 isoC'n 'PC 40:6 IsoB'n 'PI 18:0_22:5'n 'TG 50:2'n 'TG52:6'n	'TG52:3'n);*/
/****************************************************************************************************************************/
/******************************************** Exmample analysis of one figure ***********************************************/
/****************************************************************************************************************************/
libname myexcel xlsx "/home/liux3/Manuscript/26/Revised Manuscript/Data/Separate classes of data.xlsx";

data fig_1;
	set myexcel.'Fig3 A-H revised'n;

	if diet=' ' then
		delete;
run;

proc princomp data=fig_1 out=fig_1pc outstat=stats;
	* noprint;
	var 'Plasma FA 22:5'n 'Plasma TG 18:1_20:4_22:6'n 
		'Plasma TG 16:0_18:1_18:1 isomer'n 'Plasma KDdiA-PC'n 
		'Plasma Cer d 18:1_24:1'n 'Plasma Cer.d.18:1-22:0'n 
		'Plasma Acylcarnitine 16:1'n 'Plasma Cholesterol'n;
run;

data eigenvalues;
	set stats;
	where _type_='EIGENVAL';
	drop _NAME_;
run;

proc transpose data=eigenvalues out=eigen_long(rename=(COL1=eigenvalue));
	by _type_;
	var 'Plasma FA 22:5'n 'Plasma TG 18:1_20:4_22:6'n 
		'Plasma TG 16:0_18:1_18:1 isomer'n 'Plasma KDdiA-PC'n 
		'Plasma Cer d 18:1_24:1'n 'Plasma Cer.d.18:1-22:0'n 
		'Plasma Acylcarnitine 16:1'n 'Plasma Cholesterol'n;
run;

*** Apply Cheverud-Nyholt method;

proc summary data=eigen_long;
	var Eigenvalue;
	output out=stats n=M var=Var_ev;
run;

data nyholt_result;
	set stats;
	*** Cheverud-Nyholt Equation;
	M_eff2=1 + (M - 1) * (1 - (Var_ev / M));
	Alpha_adj2=0.05 / M_eff2;
	fig_id="Figure 1";
	keep fig_id M_eff2 Alpha_adj2;
run;

*** Apply Li and Ji method;

data alpha_correction;
	set eigen_long end=last;
	retain fig_id M M_eff top_pc Alpha Alpha_adj;
	retain sum_excess 0 M 0 top_pc 0;
	*** Count total variables;
	M + 1;
	*** Sum up the variance excess for all eigenvalues greater than 1;

	if eigenvalue > 1 then
		do;
			sum_excess=sum_excess + (eigenvalue - 1);
			top_pc+1;
		end;

	if last then
		do;
			M_eff=M - sum_excess;
			Alpha=0.05;
			*** Apply the PCA-based Bonferroni correction;
			Alpha_adj=Alpha / M_eff;
			*** Add a new variable here;
			fig_id="Figure 1";
			output;
		end;
	keep fig_id M M_eff Alpha Alpha_adj top_pc;
run;

data alph_correction;
	merge alpha_correction nyholt_result;
	by fig_id;

proc print;
run;

*** MANOVA analysis;
ods trace on;
ods select none;
ods output multstat=manova_p1;

proc glm data=fig_1;
	class 'germ status'n(ref='MPF') diet(ref='Palm');
	* Your categorical independent variable;
	model 'Plasma FA 22:5'n 'Plasma TG 18:1_20:4_22:6'n 
		'Plasma TG 16:0_18:1_18:1 isomer'n 'Plasma KDdiA-PC'n 
		'Plasma Cer d 18:1_24:1'n 'Plasma Cer.d.18:1-22:0'n 
		'Plasma Acylcarnitine 16:1'n 'Plasma Cholesterol'n='germ status'n diet 
		'germ status'n * diet;
	manova h=_all_/printe printh;
	* Request multivariate tests for all effects;
	run;
quit;

ods output multstat=manova_p2;

proc glm data=fig_1;
	class 'One way group'n(ref='MPF Palm');
	* Your categorical independent variable;
	model 'Plasma FA 22:5'n 'Plasma TG 18:1_20:4_22:6'n 
		'Plasma TG 16:0_18:1_18:1 isomer'n 'Plasma KDdiA-PC'n 
		'Plasma Cer d 18:1_24:1'n 'Plasma Cer.d.18:1-22:0'n 
		'Plasma Acylcarnitine 16:1'n 'Plasma Cholesterol'n='One way group'n;
	manova h=_all_/printe printh;
	* Request multivariate tests for all effects;
	run;
quit;

ods select all;
ods trace off;

data manova_p;
	set manova_p1 manova_p2;

	if hypothesis="ONE WAY GROUP" then
		hypothesis="OVERALL";
run;

data exact_p(rename=(hypothesis=factor));
	set manova_p;

	if statistic^="Wilks' Lambda" then
		delete;
	pvalue=ProbF;
	format pvalue best16.;
	keep hypothesis pvalue;
run;

proc print data=manova_p;
proc print data=exact_p;
run;
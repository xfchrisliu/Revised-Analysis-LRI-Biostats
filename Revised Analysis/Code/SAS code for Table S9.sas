options ls=200 ps=60 nodate nonumber msglevel=I validvarname=any 
	validmemname=extend;
libname mb "/home/liux3/Manuscript/26/Revised Manuscript";
libname myexcel xlsx 
	"/home/liux3/Manuscript/26/Revised Manuscript/Data/Zero 30% more.xlsx";


****** Run the Multi-Factor Monte Carlo Permutation Loop;
%macro permute_twoway(dat=, loops=);
    
    ****** 1. Clear out old results ******;
    proc datasets library=work nolist;
        delete all_perm_f;
    quit;

    ****** 2. Generate aLL shuffles upfront in one optimized DATA step ******;
    data all_shuffles;
        set &dat(keep=&var);
        
        ****** Initialize the modern random number generator ******;
        call streaminit(42);
        
        do iteration = 1 to &loops;
            rand = rand('Uniform'); * Generate a true random number *;
            output;
        end;
    run;

    ****** 3. Sort the shuffles so they are randomized within each iteration ******;
    proc sort data=all_shuffles;
        by iteration rand;
    run;

    ****** 4. Build a repetitive structural factor grid to merge against ******;
    ods select none;
    
    proc surveyselect data=&dat(keep='germ status'n Diet) 
        out=structural_grid(rename=(Replicate=iteration)) 
        method=srs 
        samprate=1 
        reps=&loops 
        seed=1234;
    run;
    
    ods select all;
    
    /*data structural_grid;
        set fig_1(keep='germ status'n Diet);
        do iteration = 1 to &loops;
            output;
        end;
    run;
    
    proc sort data=structural_grid;
        by iteration;
    run;*/

    ****** 5. Merge everything into one master simulation table ******;
    data master_iteration;
        merge structural_grid all_shuffles(keep=&var iteration);
        by iteration;
    run;

    *** Clean up massive transient files to save WORK space ***;
    proc delete data=&dat all_shuffles structural_grid; 
    quit;

    **************************************************************************** 
      6. CRITICAL OPTIMIZATION: Close ODS destinations not working in do loop
      Need to use 'ODS EXCLUDE ALL/NONE' instead of 'ODS _ALL_ CLOSE'
    ****************************************************************************;
*    ods _all_ close; 
    ods select none;

    ****** 7. Run PROC GLM exactly ONCE for all shuffles using BY processing ******;
    proc glm data=master_iteration;
        by iteration;
		class 'germ status'n(ref="MPF") diet(ref="Palm");
		model &var='germ status'n diet 'germ status'n * diet;
		ods output ModelANOVA=all_perm_f(where=(HypothesisType^=1));
    run;
    quit;

    ****** 8. Re-open ODS output so your master script can print results later ******;
*    ods html5 body="%sysfunc(pathname(work))/my_results.htm"; *%put %sysfunc(pathname(work));
    ods select all;

    ****** 9. Format outputs to exactly match what your downstream code expects ******;
    data all_perm_f;
        set all_perm_f;
        keep iteration Source FValue;
    run;

    *** Clean up the large master table ***;
    proc delete data=master_iteration;
    quit;

%mend permute_twoway;


****** Alternative way to generate shuffled data and results ******;
/*%macro permute_twoway(dat=, loops=);
	
	****** Clear out any old appended datasets ********;
	proc datasets library=work nolist;
		delete all_perm_f;
		run;
	quit;

	%do k=1 %to &loops;
		ods select none;
		ods graphics off;

		****** Shuffle the raw variable values randomly ******;
		data shuffled;
			set &dat(keep=&var);
			rand=ranuni(42 + &k);
		run;

		proc sort data=shuffled;
			by rand;
		run;

		****** Merge shuffled values back onto fixed structural factors ******;
		data iteration;
			merge &dat(keep='germ status'n Diet) shuffled(keep=&var);
		run;

		****** Calculate F-stats for this specific shuffle ******;
		proc glm data=iteration;
			class 'germ status'n(ref="MPF") diet(ref="Palm");
			model &var='germ status'n diet 'germ status'n * diet;
			ods output ModelANOVA=f_out(where=(HypothesisType^=1));
			run;

		data f_out;
			set f_out;
			iteration=&k;
			keep iteration Source FValue;
		run;

		proc append base=all_perm_f data=f_out force;
			run;;

		proc delete data=&dat shuffled iteration f_out;
		run;

	%end;
	ods graphics on;
	ods select all;
%mend permute_twoway;*/

****** For one variable in a figure;
%macro one_variable (fig, var, nperm);
	data fig_1;
		set myexcel.&fig;
		if diet=' ' then delete;
	run;

	****** Calculate observed F-statistics for Main Effects and Interaction ******;
	ods select none;

	proc glm data=fig_1;
		class 'germ status'n(ref="MPF") diet(ref="Palm");
		model &var='germ status'n diet 'germ status'n * diet;
		ods output ModelANOVA=observed_f;
		run;
		
	ods select all;

    ******* Filter down to just the effects we care about ******;
	data observed_f(keep=Source Obs_F);
		set observed_f;
		format Fvalue best10.;
		if HypothesisType=1 then delete;
		rename Fvalue=Obs_F;
	run;

	****** Run the Multi-Factor Monte Carlo Permutation Loop ******;
	%permute_twoway(dat=fig_1, loops=&nperm);

	proc sort data=observed_f;
		by Source;
	run;

	proc sort data=all_perm_f;
		by Source;
	run;

	/*data final_two_way_p;
	merge all_perm_f observed_f;
	by Source;
	run;*/
	data final_p;
		merge all_perm_f observed_f;
		by Source;
		retain greater_or_equal 0 total_shuffles 0;

		if first.Source then
			do;
				greater_or_equal=0;
				total_shuffles=0;
			end;

		total_shuffles + 1;

		if FValue >=Obs_F then
			greater_or_equal + 1;

		if last.Source then
			do;
				Permutation_P=greater_or_equal / total_shuffles;
				output;
			end;
		*    keep Source Obs_F Permutation_P;
	run;

	proc transpose data=final_p out=wide_final_p;
		id source;
		var Permutation_P;
	run;

	data pvalue_AllZero;
		retain fig_id variable 'germ status'n diet;
		set wide_final_p;
		length fig_id $25 variable $30;
		fig_id="%sysfunc(dequote(&fig))";
		;
		variable="%sysfunc(dequote(&var))";
		drop _name_;
		format 'germ status'n diet 'germ status*diet'n best10.;
	run;

	/*proc print data=pvalue_AllZero;
	run;*/

%mend one_variable;

/*%one_variable('Fig 4 A-H'n, 'LIVER aKETOD'n, 2000);
%one_variable('Fig 4 I-P'n, 'LIVER RESOLVIN D5'n, 2000);
%one_variable('Fig S4 A-H'n, 'PLASMA 12-OXO'n, 2000);
%one_variable('Fig S4 I-P'n, 'ILEUM aKETOD'n, 2000);
%one_variable('Fig S4 I-P'n, 'ILEUM gHYD'n, 2000);*/


****** For all variables and figures ******;
%macro all_figure(fad, nsample);
	%local num_fad a this_ad num_fid i this_id fig;
	%let num_fad=%sysfunc(countw(&fad, %str(~)));

	%do a=1 %to &num_fad;
		%let this_ad=%scan(&fad, &a, %str(~));
		%let num_fid=%sysfunc(countw(&this_ad, %str(|)));

		%do i=1 %to &num_fid;
			%let this_id=%scan(&this_ad, &i, %str(|));

			%if &i=1 %then
				%let fig=&this_id;
			%else
				%do;
					%one_variable(&fig, &this_id, &nsample);

					%if &i=2 %then
						%do;

							data pvalue_AllZero&a;
								set pvalue_AllZero;
							run;

						%end;
					%else
						%do;

							data pvalue_AllZero&a;
								set pvalue_AllZero&a pvalue_AllZero;
							run;

						%end;
				%end;
		%end;
	%end;

	data pvalue_AllZero_all;
		set pvalue_AllZero1-pvalue_AllZero&num_fad;
	proc print;
	run;

%mend all_figure;


***** 30% more zeros;
%all_figure('Fig 4 A-H'n|'LIVER aHYA'n|'LIVER aKETOD'n|'LIVER HYE'n|'LIVER gHYA'n~
'Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n|'LIVER RESOLVIN D5'n|'LIVER 20-HETE'n|'LIVER N-DGLA EA'n~
'Fig S5 A-H revised done'n|'PLASMA aHYA+aHYC'n|'PLASMA gHYD'n|'PLASMA 12-HAS'n|'PLASMA 12-OXO'n~
'Fig S5 I-P revised done'n|'ILEUM aKETOD'n|'ILEUM gHYD'n ~
'Fig S6 revised done'n|'Leukotriene B4'n|'15-deoxy-PGJ2'n|'N-C22:6n3_EA'n|'13-KODE'n~, 10000);

proc export data=pvalue_AllZero_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_Zero30_all.csv" 
		dbms=csv replace;
run;


****** 50% more zeros;
%all_figure('Fig 4 A-H'n|'LIVER aHYA'n|'LIVER aKETOD'n|'LIVER gHYA'n~
'Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n|'LIVER RESOLVIN D5'n|'LIVER 20-HETE'n|'LIVER N-DGLA EA'n~
'Fig S5 A-H revised done'n|'PLASMA aHYA+aHYC'n|'PLASMA 12-HAS'n|'PLASMA 12-OXO'n~
'Fig S5 I-P revised done'n|'ILEUM aKETOD'n|'ILEUM gHYD'n ~
'Fig S6 revised done'n|'Leukotriene B4'n|'15-deoxy-PGJ2'n~, 10000);

proc export data=pvalue_AllZero_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_Zero50_all.csv" 
		dbms=csv replace;
run;


****** 60% more zeros;
%all_figure('Fig 4 A-H'n|'LIVER aHYA'n|'LIVER aKETOD'n|'LIVER gHYA'n ~
'Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n| 'LIVER RESOLVIN D5'n ~
'Fig S5 A-H revised done'n|'PLASMA aHYA+aHYC'n| 'PLASMA 12-OXO'n ~
'Fig S5 I-P revised done'n|'ILEUM aKETOD'n~
'Fig S6 revised done'n|'Leukotriene B4'n|'15-deoxy-PGJ2'n ~, 100);

proc export data=pvalue_AllZero_all 
		outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_Zero60_all.csv" 
		dbms=csv replace;
run;





/*************************************************************************************************************************/
/********************* Example of permutation test for two factors and one variable/metabolite ***************************/
/*************************************************************************************************************************/
libname myexcel xlsx 
	"/home/liux3/Manuscript/26/Mark Brown/Data/Updated Data/Zero 30% more.xlsx";

****** Fig S4 A-H and variable PLASMA 12-OXO;
data fig_1;
	set myexcel.'Fig S4 A-H'n;
	if diet=' ' then delete;
run;

****** 1. Calculate observed F-statistics for Main Effects and Interaction ******;
ods select none;

proc glm data=fig_1;
	class 'germ status'n(ref="MPF") diet(ref="Palm");
	model 'PLASMA 12-OXO'n='germ status'n diet 'germ status'n * diet;
	ods output ModelANOVA=observed_f;
run;

ods select all;

****** Filter down to just the effects we care about ******;
data observed_f(keep=Source Obs_F);
	set observed_f;
	if HypothesisType=1 then
		delete;
	rename Fvalue=Obs_F;
run;

****** 2. Run the Multi-Factor Monte Carlo Permutation Loop ******;
%macro permute_twoway(var, dat=, loops=);

	*** Clear out any old appended datasets ***;
	proc datasets library=work nolist;
		delete all_perm_f;
		run;
	quit;

	%do i=1 %to &loops;
		ods select none;

		*** Shuffle the raw variable values randomly ***;
		data shuffled;
			set &dat(keep=&var);
			rand=ranuni(42 + &i);
		run;

		proc sort data=shuffled;
			by rand;
		run;

		*** Merge shuffled values back onto fixed structural factors ***;
		data iteration;
			merge &dat(keep='germ status'n Diet) shuffled(keep=&var);
		run;

		*** Calculate F-stats for this specific shuffle ***;
		proc glm data=iteration;
			class 'germ status'n(ref="MPF") diet(ref="Palm");
			model &var='germ status'n diet 'germ status'n * diet;
			ods output ModelANOVA=f_out(where=(HypothesisType^=1));
			run;

		data f_out;
			set f_out;
			iteration=&i;
			keep iteration Source FValue;
		run;

		proc append base=all_perm_f data=f_out force;
		run;

		proc delete data=&dat shuffled iteration f_out;
		run;

	%end;
	ods select all;
%mend permute_twoway;

%permute_twoway('PLASMA 12-OXO'n, dat=fig_1, loops=2000);

proc sort data=observed_f;
	by Source;
run;

proc sort data=all_perm_f;
	by Source;
run;

/*data final_two_way_p;
merge all_perm_f observed_f;
by Source;
run;*/
data final_p;
	merge all_perm_f observed_f;
	by Source;
	retain greater_or_equal 0 total_shuffles 0;

	if first.Source then
		do;
			greater_or_equal=0;
			total_shuffles=0;
		end;
	total_shuffles + 1;

	if FValue >=Obs_F then
		greater_or_equal + 1;

	if last.Source then
		do;
			Permutation_P=greater_or_equal / total_shuffles;
			output;
		end;
	*    keep Source Obs_F Permutation_P;
run;

proc transpose data=final_p out=wide_final_p;
	id source;
	var Permutation_P;
run;

data pvalue_AllZero;
	retain fig_id variable 'GERM STATUS'n DIET;
	set wide_final_p;
	length fig_id $25 variable $30;
	fig_id="Fig S4 A-H";
	variable="PLASMA 12-OXO";
	drop _name_;
	format 'germ status'n diet 'germ status*diet'n best10.;
run;

*********************************************************************************************************;
*********************************************************************************************************;
libname myexcel xlsx 
	"/home/liux3/Manuscript/26/Mark Brown/Data/Updated Data/AllZero_Group 2.xlsx";

****** Fig 4 A-H and variable LIVER aKETOD;
data fig_1;
	set myexcel.'Fig 4 A-H'n;
	if diet=' ' then delete;
run;

/* 1. Calculate observed F-statistics for Main Effects and Interaction */
ods select none;

proc glm data=fig_1;
	class 'germ status'n(ref="MPF") diet(ref="Palm");
	model 'LIVER aKETOD'n='germ status'n diet 'germ status'n * diet;
	ods output ModelANOVA=observed_f;
	run;
	ods select all;

	/* Filter down to just the effects we care about */
data observed_f(keep=Source Obs_F);
	set observed_f;
	if HypothesisType=1 then delete;
	rename Fvalue=Obs_F;
run;

/* 2. Run the Multi-Factor Monte Carlo Permutation Loop */
%macro permute_twoway(var, dat=, loops=);

	/* Clear out any old appended datasets */
	proc datasets library=work nolist;
		delete all_perm_f;
	run;
	quit;

	%do i=1 %to &loops;
		ods select none;

		/* Shuffle the raw variable values randomly */
		data shuffled;
			set &dat(keep=&var);
			rand=ranuni(42 + &i);
		run;

		proc sort data=shuffled;
			by rand;
		run;

		/* Merge shuffled values back onto fixed structural factors */
		data iteration;
			merge &dat(keep='germ status'n Diet) shuffled(keep=&var);
		run;

		/* Calculate F-stats for this specific shuffle */
		proc glm data=iteration;
			class 'germ status'n(ref="MPF") diet(ref="Palm");
			model &var='germ status'n diet 'germ status'n * diet;
			ods output ModelANOVA=f_out(where=(HypothesisType^=1));
			run;

		data f_out;
			set f_out;
			iteration=&i;
			keep iteration Source FValue;
		run;

		proc append base=all_perm_f data=f_out force;
		run;

		proc delete data=&dat shuffled iteration f_out;
		run;

	%end;
	ods select all;
%mend permute_twoway;

%permute_twoway('LIVER aKETOD'n, data=fig_1, loops=2000);

proc sort data=observed_f;
	by Source;
run;

proc sort data=all_perm_f;
	by Source;
run;


data final_p;
	merge all_perm_f observed_f;
	by Source;
	retain greater_or_equal 0 total_shuffles 0;

	if first.Source then
		do;
			greater_or_equal=0;
			total_shuffles=0;
		end;
	total_shuffles + 1;

	if FValue >=Obs_F then
		greater_or_equal + 1;

	if last.Source then
		do;
			Permutation_P=greater_or_equal / total_shuffles;
			output;
		end;
	*    keep Source Obs_F Permutation_P;
run;

proc transpose data=final_p out=wide_final_p;
	id source;
	var Permutation_P;
run;

data pvalue_AllZero;
	retain fig_id variable 'germ status'n diet;
	set wide_final_p;
	length fig_id $25 variable $30;
	fig_id="Fig 4 A-H";
	variable="LIVER aKETOD";
	drop _name_;
	format 'germ status'n diet 'germ status*diet'n best10.;
run;

/*********************** Permutatin test for one factor with three groups ***************************/
/*data diet_data;
input Group $ Value;
datalines;
Diet1 0
Diet1 0
Diet1 0
Diet2 0
Diet2 0
Diet2 0
Diet3 5
Diet3 7
Diet3 6
;

proc npar1way data=diet_data SCORES=DATA;
class Group;
var Value;
exact scores=data / mc n=10000 seed=42;
run;*/
/*proc multtest data=diet_data permutation nsample=10000 seed=1234;
class Group;
test mean(Value);
run;*/
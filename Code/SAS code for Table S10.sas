options ls=200 ps=60 nodate nonumber msglevel=I validvarname=any 
	validmemname=extend;
libname xf "/home/liux3/Manuscript/26/Revised Manuscript";
libname myexcel xlsx 
	"/home/liux3/Manuscript/26/Revised Manuscript/Data/Zero 30% more.xlsx";

/* =====================================================================
1. TWO-WAY MONTE CARLO PERMUTATION LOOP
===================================================================== */
%macro permute_twoway(dat=, loops=);
    
    ****** Clear out old results ******;
    proc datasets library=work nolist;
        delete all_perm_germ; *Diet all_perm_diet all_perm_germdiet;
    quit;
    
    options nonotes nosource nomprint;
    ods select none;
    ods graphics off;

    ****** Generate ALL shuffleed datasets upfront in memory ******;
    data all_shuffles;
        set &dat(keep=&var);
        
        ****** Initialize the modern random number generator ******;
        call streaminit(42);      
        do iteration = 1 to &loops;
            rand = rand('Uniform'); * Generate a true random number *;
            output;
        end;
    run;
    
    ****** Sort the shuffles so they are randomized within each iteration ******;
    proc sort data=all_shuffles;
        by iteration rand;
    run;

    ****** Build a repetitive structural factor grid to merge against ******;    
    proc surveyselect data=&dat(keep='germ status'n Diet) 
        out=structural_grid(rename=(Replicate=iteration)) 
        method=srs 
        samprate=1 
        reps=&loops 
        seed=1234;
    run;
      
    ****** Merge everything into one master simulation table ******;
    data master_iteration;
        merge structural_grid all_shuffles(keep=&var iteration);
        by iteration;
    run;

    *** Clean up massive transient files to save WORK space ***;
    proc delete data=&dat all_shuffles structural_grid; 
    quit;

    ****** Run PROC GLM exactly ONCE for all shuffles using BY processing ******;
    proc glm data=master_iteration;
        by iteration;
        class 'germ status'n(ref="MPF") diet(ref="Palm");
        model &var='germ status'n diet 'germ status'n * diet;
        lsmeans 'germ status'n*diet /slice=diet tdiff pdiff=control('MPF' 'Palm') adjust=dunnett;
        
        *** Direct output routing bypassing ODS memory layers ***;
        ods output SlicedANOVA=all_perm_germ; *LSMeans=all_perm_germdiet SlicedANOVA=all_perm_germ;
    run;
    quit;

    /*proc sort data=master_iteration;
        by iteration 'germ status'n;
    run;
    
    proc glm data=master_iteration;
        by iteration 'germ status'n;
        class diet(ref="Palm");
        model &var=diet;
        lsmeans diet /tdiff pdiff=control('Palm') adjust=dunnett;
        ods output LSMeans=all_perm_diet;
    run;
    quit;*/

*    options notes source mprint;
    ods select all;
    ods graphics on;

    ****** Format outputs to exactly match what your downstream code expects ******;
    data all_perm_germ;
        set all_perm_germ;
        format FValue best10.;
        keep iteration Diet FValue;
    run;

    /*data all_perm_diet;
        set all_perm_diet;
        keep iteration 'Germ Status'n Diet tDiff;
    run;

    data all_perm_germdiet;
        set all_perm_germdiet;
        keep iteration 'Germ Status'n Diet tDiff;
    run;*/

    *** Clean up the large master table ***;
    proc delete data=master_iteration;
    quit;

%mend permute_twoway;


******************************************************************** 
       Alternative way to generate shuffled data and results. 
       This is an algorithm that may work well as expected 
********************************************************************;
/*%macro permute_twoway(dat=, loops=);
    
    ****** 1. Clear out old results ******;
    proc datasets library=work nolist;
        delete all_perm_germ all_perm_diet all_perm_germdiet;
    quit;

    ****** 2. Generate ALL shuffles upfront in one optimized DATA step ******;
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
    data structural_grid;
        set &dat(keep='germ status'n Diet);
        do iteration = 1 to &loops;
            output;
        end;
    run;
    
    proc sort data=structural_grid;
        by iteration;
    run;

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
    ****************************************************************************;
*    ods _all_ close; 
    ods select none;

    ****** 7. Run PROC GLM exactly ONCE for all shuffles using BY processing ******;
    proc glm data=master_iteration;
        by iteration;
        class 'germ status'n(ref="MPF") diet(ref="Palm");
        model &var='germ status'n diet 'germ status'n * diet;
        lsmeans 'germ status'n*diet /slice=diet tdiff pdiff=control('MPF' 'Palm') adjust=dunnett;
        
        *** Direct output routing bypassing ODS memory layers ***;
        ods output LSMeans=all_perm_germdiet SlicedANOVA=all_perm_germ;
    run;
    quit;

    proc sort data=master_iteration;
        by iteration 'germ status'n;
    run;

    proc glm data=master_iteration;
        by iteration 'germ status'n;
        class diet(ref="Palm");
        model &var=diet;
        lsmeans diet /tdiff pdiff=control('Palm') adjust=dunnett;
        ods output LSMeans=all_perm_diet;
    run;
    quit;

    ****** 8. Re-open ODS output so your master script can print results later ******;
*    ods html5 body="%sysfunc(pathname(work))/my_results.htm"; *%put %sysfunc(pathname(work));
    ods select all;

    ****** 9. Format outputs to exactly match what your downstream code expects ******;
    data all_perm_germ;
        set all_perm_germ;
        keep iteration Diet FValue;
    run;

    data all_perm_diet;
        set all_perm_diet;
        keep iteration 'Germ Status'n Diet tDiff;
    run;

    data all_perm_germdiet;
        set all_perm_germdiet;
        keep iteration 'Germ Status'n Diet tDiff;
    run;

    *** Clean up the large master table ***;
    proc delete data=master_iteration;
    quit;

%mend permute_twoway;*/


***********************************************************************
       Alternative way to generate shuffled data and results using
       a do loop that takes much longer to get the final results
***********************************************************************;
/*%macro permute_twoway(dat=, loops=);
	
	****** Clear out any old appended datasets ******;
	proc datasets library=work nolist;
		delete all_perm_germ all_perm_diet all_perm_germdiet;
		run;
	quit;

	%do k=1 %to &loops;
		ods select none;

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

		****** Calculate F or t-stats for this specific shuffle ******;
		proc glm data=iteration;
			class 'germ status'n(ref="MPF") diet(ref="Palm");
			model &var='germ status'n diet 'germ status'n * diet;
			lsmeans 'germ status'n*diet /slice=diet tdiff pdiff=control('MPF' 'Palm') 
				adjust=dunnett;
			ods output LSMeans=germdiett_out SlicedANOVA=germf_out;
			run;

		proc sort data=iteration;
			by 'germ status'n;
		run;

		proc glm data=iteration;
			by 'germ status'n;
			class diet(ref="Palm");
			model &var=diet;
			lsmeans diet /tdiff pdiff=control('Palm') adjust=dunnett;
			ods output LSMeans=diett_out;
			run;

		data germf_out;
			set germf_out;
			iteration=&k;
			keep iteration Diet FValue;
		run;

		data diett_out;
			set diett_out;
			iteration=&k;
			keep iteration 'Germ Status'n Diet tDiff;
		run;

		data germdiett_out;
			set germdiett_out;
			iteration=&k;
			keep iteration 'Germ Status'n Diet tDiff;
		run;

		****** Combined all the simulated output results ******;
		proc append base=all_perm_germ data=germf_out force;
		run;

		proc append base=all_perm_diet data=diett_out force;
		run;

		proc append base=all_perm_germdiet data=germdiett_out force;
		run;

		proc delete data=shuffled iteration germf_out diett_out germdiett_out;
		run;
		
	%end;
	
	ods select all;
	
%mend permute_twoway;*/

/* =====================================================================
         2. PROCESS ONE VARIABLE WITHIN A FIGURE
   ===================================================================== */
%macro one_variable(fig, var, nperm);
	data fig_1;
		set myexcel.&fig;
		if diet=' ' then delete;
	run;

	****** Calculate observed F-statistics for Main Effects and Interaction ******;
	ods select none;

	proc glm data=fig_1;
		class 'germ status'n(ref="MPF") diet(ref="Palm");
		model &var='germ status'n diet 'germ status'n * diet;
		lsmeans 'germ status'n*diet /slice=diet tdiff pdiff=control('MPF' 'Palm') 
			adjust=dunnett;
		ods output SlicedANOVA=germdiff; *LSMeans=germdietdiff SlicedANOVA=germdiff;
		run;
	quit;

    /*proc sort data=fig_1;
		by 'germ status'n;
	run;

	proc glm data=fig_1;
		by 'germ status'n;
		class diet(ref="Palm");
		model &var=diet;
		lsmeans diet /tdiff pdiff=control('Palm') adjust=dunnett;
		ods output LSMeans=dietdiff;
		run;
	quit;*/

	****** Filter down to just the effects we care about ******;
	data germdiff (keep=Diet germObs_F);
		set germdiff;
		rename Fvalue=germObs_F;
	run;

    /*data dietdiff (keep='Germ Status'n Diet dietObs_t);
		set dietdiff;
		rename tDiff=dietObs_t;
	run;

	data germdietdiff (keep='Germ Status'n Diet germdietObs_t);
		set germdietdiff;
		rename tDiff=germdietObs_t;
	run;*/

	ods select all;

	****** Run the Multi-Factor Monte Carlo Permutation Loop ******;
    %permute_twoway(dat=fig_1, loops=&nperm);


	****** p values for germ given each level of diet ******;
	proc sort data=germdiff;
		by Diet;
	run;

	proc sort data=all_perm_germ;
		by Diet;
	run;

	data final_pgerm;
		merge all_perm_germ germdiff;
		by Diet;
		retain greater_or_equal 0 total_shuffles 0;

		if first.Diet then
			do;
				greater_or_equal=0;
				total_shuffles=0;
			end;
		total_shuffles + 1;

		if FValue >=germObs_F then
			greater_or_equal + 1;

		if last.Diet then
			do;
				Permutation_P=greater_or_equal / total_shuffles;
				output;
			end;
		*    keep Source Obs_F Permutation_P;
	run;

	proc transpose data=final_pgerm out=wide_final_pgerm;
		id Diet;
		var Permutation_P;
	run;

	data pvalue_AllZero_germ;
		retain fig_id variable;
		set wide_final_pgerm;
		length fig_id $25 variable $30;
		fig_id="%sysfunc(dequote(&fig))";
		variable="%sysfunc(dequote(&var))";
		drop _name_;
		format Borage Echium Fish HLS HOS Palm best10.;
	run;
	
	proc delete data=germdiff all_perm_germ final_pgerm wide_final_pgerm;
	quit;


	****** p values for diet given each level of germ status ******;
    /*proc sort data=dietdiff;
		by 'Germ Status'n Diet;
	run;

	proc sort data=all_perm_diet;
		by 'Germ Status'n Diet;
	run;

	data final_pdiet;
		merge all_perm_diet dietdiff;
		by 'Germ Status'n Diet;
		retain greater_or_equal 0 total_shuffles 0;

		if first.Diet then
			do;
				greater_or_equal=0;
				total_shuffles=0;
			end;
		total_shuffles + 1;

		if abs(tDiff) >=abs(dietObs_t) then
			greater_or_equal + 1;

		if last.Diet then
			do;
				Permutation_P=greater_or_equal / total_shuffles;
				output;
			end;
		*    keep Source Obs_F Permutation_P;
	run;

	proc transpose data=final_pdiet out=wide_final_pdiet;
		by 'Germ Status'n;
		id Diet;
		var Permutation_P;
	run;

	data pvalue_AllZero_diet;
		retain fig_id variable;
		set wide_final_pdiet;
		length fig_id $25 variable $30;
		fig_id="%sysfunc(dequote(&fig))";
		variable="%sysfunc(dequote(&var))";
		drop _name_ Palm;
		format Borage Echium Fish HLS HOS best10.;
	run;
	
	proc delete data=dietdiff all_perm_diet final_pdiet wide_final_pdiet;
	quit;


	****** p values for germ status and diet combinations ******;
	proc sort data=germdietdiff;
		by 'Germ Status'n Diet;
	run;

	proc sort data=all_perm_germdiet;
		by 'Germ Status'n Diet;
	run;

	data final_pgermdiet;
		merge all_perm_germdiet germdietdiff;
		by 'Germ Status'n Diet;
		retain greater_or_equal 0 total_shuffles 0;

		if first.Diet then
			do;
				greater_or_equal=0;
				total_shuffles=0;
			end;
		total_shuffles + 1;

		if abs(tDiff) >=abs(germdietObs_t) then
			greater_or_equal + 1;

		if last.Diet then
			do;
				Permutation_P=greater_or_equal / total_shuffles;
				output;
			end;
		*    keep Source Obs_F Permutation_P;
	run;

	proc transpose data=final_pgermdiet out=wide_final_pgermdiet;
		id 'Germ Status'n Diet;
		var Permutation_P;
	run;

	data pvalue_AllZero_germdiet;
		retain fig_id variable;
		set wide_final_pgermdiet;
		length fig_id $25 variable $30;
		fig_id="%sysfunc(dequote(&fig))";
		variable="%sysfunc(dequote(&var))";
		drop _name_ MPFPalm;
		format GFBorage GFEchium GFFish GFHLS GFHOS GFPalm MPFBorage MPFEchium MPFFish MPFHLS MPFHOS best10.;
	run;

	proc delete data=germdietdiff all_perm_germdiet final_pgermdiet wide_final_pgermdiet;
	quit;*/

%mend one_variable;

/* =====================================================================
         3. MASTER MACRO TO EXECUTE ALL FIGURES & VARIABLES
   ===================================================================== */
%macro all_figure(fad, nsample);
	%local num_fad a this_ad num_fid i this_id fig;
    ods exclude all;
    
	****** Clean out master long file from any old sessions ******;
	proc datasets library=work nolist;
		delete pvalue_AllZero_germ_all; *pvalue_AllZero_diet_all pvalue_AllZero_germdiet_all;
	run;
	quit;

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

				 ****** Call the single variable processor ******;
                 %one_variable(&fig, &this_id, &nsample);

					****** Append wide format variables immediately to a master dataset ******;
					proc append base=pvalue_AllZero_germ_all data=pvalue_AllZero_germ force;
					run;

					/*proc append base=pvalue_AllZero_diet_all data=pvalue_AllZero_diet force;
					run;

					proc append base=pvalue_AllZero_germdiet_all data=pvalue_AllZero_germdiet 
							force;
					run;*/

					proc delete data=pvalue_AllZero_germ; *pvalue_AllZero_diet 
							pvalue_AllZero_germdiet;
					run;

				%end;
		%end;
	%end;
	
	ods exclude none;

	****** Print out the consolidated clean results matrix ******;
	proc print data=pvalue_AllZero_germ_all;
		title "Consolidated Permutation P-Values for Germ";
	run;

	title;

    /*proc print data=pvalue_AllZero_diet_all noobs;
		title "Consolidated Permutation P-Values for Diet";
	run;

	title;

	proc print data=pvalue_AllZero_germdiet_all noobs;
		title "Consolidated Permutation P-Values for Germ and diet";
	run;

	title;*/
%mend all_figure;

/* =====================================================================
4. EXECUTION
===================================================================== */

***** 30% more zeros;
%all_figure('Fig 4 A-H'n|'LIVER aHYA'n|'LIVER aKETOD'n|'LIVER HYE'n|'LIVER gHYA'n~
'Fig 4 I-P'n|'LIVER 18-RESOLVIN E2'n|'LIVER RESOLVIN D5'n|'LIVER 20-HETE'n|'LIVER N-DGLA EA'n~
'Fig S5 A-H revised done'n|'PLASMA aHYA+aHYC'n|'PLASMA gHYD'n|'PLASMA 12-HAS'n|'PLASMA 12-OXO'n~
'Fig S5 I-P revised done'n|'ILEUM aKETOD'n|'ILEUM gHYD'n ~
'Fig S6 revised done'n|'Leukotriene B4'n|'15-deoxy-PGJ2'n|'N-C22:6n3_EA'n|'13-KODE'n~, 10000);

proc export data=pvalue_AllZero_germ_all 
            outfile="/home/liux3/Manuscript/26/Revised Manuscript/Results/pvalue_Zero30_germ_all.csv" 
            dbms=csv 
            replace;
run;


/*proc export data=pvalue_AllZero_diet_all 
            outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_Zero30_diet_all.csv" 
            dbms=csv 
            replace;
run;

proc export data=pvalue_AllZero_germdiet_all 
            outfile="/home/liux3/Manuscript/26/Mark Brown/Results/pvalue_Zero30_germdiet_all.csv" 
            dbms=csv 
            replace;
run;*/





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

****** Calculate observed F- or t-statistics for Simple Effects ******;
ods select none;

proc glm data=fig_1;
	class 'germ status'n(ref="MPF") diet(ref="Palm");
	model 'PLASMA 12-OXO'n='germ status'n diet 'germ status'n * diet;
	lsmeans 'germ status'n*diet /slice=diet tdiff pdiff=control('MPF' 'Palm') 
		adjust=dunnett;
	ods output SlicedANOVA=germdiff; *LSMeans=germdietdiff SlicedANOVA=germdiff;
	run;
quit;

****** Filter down to just the effects we care about ******;
data germdiff (keep=Diet germObs_F);
	set germdiff;
	format Fvalue best10.;
	rename Fvalue=germObs_F;
run;

ods select all;

****** Run the Multi-Factor Monte Carlo Permutation Loop ******;
%macro permute_twoway(var, dat=, loops=);
    
    ****** Clear out old results ******;
    proc datasets library=work nolist;
        delete all_perm_germ; *Diet all_perm_diet all_perm_germdiet;
    quit;
    
    options nonotes nosource nomprint;
    ods select none;
    ods graphics off;

    ****** Generate ALL shuffleed datasets upfront in memory ******;
    data all_shuffles;
        set &dat(keep=&var);
        
        ****** Initialize the modern random number generator ******;
        call streaminit(42);       
        do iteration = 1 to &loops;
            rand = rand('Uniform'); * Generate a true random number *;
            output;
        end;
    run;
    
    ****** Sort the shuffles so they are randomized within each iteration ******;
    proc sort data=all_shuffles;
        by iteration rand;
    run;

    ****** Build a repetitive structural factor grid to merge against ******;    
    proc surveyselect data=&dat(keep='germ status'n Diet) 
        out=structural_grid(rename=(Replicate=iteration)) 
        method=srs 
        samprate=1 
        reps=&loops 
        seed=1234;
    run;
      
    ****** Merge everything into one master simulation table ******;
    data master_iteration;
        merge structural_grid all_shuffles(keep=&var iteration);
        by iteration;
    run;

    *** Clean up massive transient files to save WORK space ***;
    proc delete data=&dat all_shuffles structural_grid; 
    quit;

    ****** Run PROC GLM exactly ONCE for all shuffles using BY processing ******;
    proc glm data=master_iteration;
        by iteration;
        class 'germ status'n(ref="MPF") diet(ref="Palm");
        model &var='germ status'n diet 'germ status'n * diet;
        lsmeans 'germ status'n*diet /slice=diet tdiff pdiff=control('MPF' 'Palm') adjust=dunnett;
        
        *** Direct output routing bypassing ODS memory layers ***;
        ods output SlicedANOVA=all_perm_germ; *Diet LSMeans=all_perm_germdiet SlicedANOVA=all_perm_germ;
    run;
    quit;

    /*proc sort data=master_iteration;
        by iteration 'germ status'n;
    run;
    
    proc glm data=master_iteration;
        by iteration 'germ status'n;
        class diet(ref="Palm");
        model &var=diet;
        lsmeans diet /tdiff pdiff=control('Palm') adjust=dunnett;
        ods output LSMeans=all_perm_diet;
    run;
    quit;*/

*    options notes source mprint;
    ods select all;
    ods graphics on;

    ****** Format outputs to exactly match what your downstream code expects ******;
    data all_perm_germ;
        set all_perm_germ;
        keep iteration Diet FValue;
    run;

    /*data all_perm_diet;
        set all_perm_diet;
        keep iteration 'Germ Status'n Diet tDiff;
    run;

    data all_perm_germdiet;
        set all_perm_germdiet;
        keep iteration 'Germ Status'n Diet tDiff;
    run;*/

    *** Clean up the large master table ***;
    proc delete data=master_iteration;
    quit;

%mend permute_twoway;

%permute_twoway('PLASMA 12-OXO'n, dat=fig_1, loops=10000);


****** Alternative way to generate shuffled data and results ******;
/*%macro permute_twoway(var, dat=, loops=);
    
    ****** 1. Clear out old results ******;
    proc datasets library=work nolist;
        delete all_perm_germ ;
    quit;

    ****** 2. Generate ALL shuffles upfront in one optimized DATA step ******;
    data all_shuffles;
        set &dat(keep=&var);
        do iteration = 1 to &loops;
            *** Generate a unique random number for every iteration ***;
            rand = ranuni(42 + iteration); 
            output;
        end;
    run;

    ****** 3. Sort the shuffles so they are randomized within each iteration ******;
    proc sort data=all_shuffles;
        by iteration rand;
    run;

    ****** 4. Build a repetitive structural factor grid to merge against ******;
    data structural_grid;
        set &dat(keep='germ status'n Diet);
        do iteration = 1 to &loops;
            output;
        end;
    run;
    
    proc sort data=structural_grid;
        by iteration;
    run;

    ****** 5. Merge everything into one master simulation table ******;
    data master_iteration;
        merge structural_grid all_shuffles(keep=&var iteration);
        by iteration;
    run;

    *** Clean up massive transient files to save WORK space ***;
    proc delete data=&dat all_shuffles structural_grid; 
    quit;

    ****** 6. CRITICAL OPTIMIZATION: Close ODS destinations entirely ******;
*    ods _all_ close; 
    ods select none;

    ****** 7. Run PROC GLM exactly ONCE for all shuffles using BY processing ******;
    proc glm data=master_iteration;
        by iteration;
        class 'germ status'n(ref="MPF") diet(ref="Palm");
        model &var='germ status'n diet 'germ status'n * diet;
        lsmeans 'germ status'n*diet /slice=diet tdiff pdiff=control('MPF' 'Palm') adjust=dunnett;
        
        *** Direct output routing bypassing ODS memory layers ***;
        ods output SlicedANOVA=all_perm_germ;
    run;
    quit;

    ****** 8. Re-open ODS output so your master script can print results later ******;
*    ods html5 body="%sysfunc(pathname(work))/my_results.htm"; 
    ods select all;

    ****** 9. Format outputs to exactly match what your downstream code expects ******;
    data all_perm_germ;
        set all_perm_germ;
        keep iteration Diet FValue;
    run;

    *** Clean up the large master table ***;
    proc delete data=master_iteration;
    quit;

%mend permute_twoway;

%permute_twoway('PLASMA 12-OXO'n, dat=fig_1, loops=100);*/


****** p values for germ given each level of diet ******;
proc sort data=germdiff;
	by Diet;
run;

proc sort data=all_perm_germ;
	by Diet;
run;

data final_pgerm;
	merge all_perm_germ germdiff;
	by Diet;
	retain greater_or_equal 0 total_shuffles 0;

	if first.Diet then
		do;
			greater_or_equal=0;
			total_shuffles=0;
		end;
	total_shuffles + 1;

	if FValue >=germObs_F then
		greater_or_equal + 1;

	if last.Diet then
		do;
			Permutation_P=greater_or_equal / total_shuffles;
			output;
		end;
	*    keep Source Obs_F Permutation_P;
run;

proc transpose data=final_pgerm out=wide_final_pgerm;
	id Diet;
	var Permutation_P;
run;

data pvalue_AllZero_germ;
	retain fig_id variable;
	set wide_final_pgerm;
	length fig_id $25 variable $30;
	fig_id="Fig S4 A-H";
	variable="PLASMA 12-OXO";
	drop _name_;
	*format ;
run;
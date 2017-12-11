/****************************************************************************

-----------------------------------------------------------------------------
--------------------           SiGe_TO_dist_angles         ------------------
----              Univ. Pablo Olavide, RASPA Group, Spain                ----
------------------------        Janauary 2017       -------------------------
-----------------------------------------------------------------------------

  utilitary program to extract TO distances, OTO angles, TOT angles and T-T
  distances for a zeolite 
 
  crash_label is an error message flag: = 0 if error ocurrs, > 0 else 
*****************************************************************************/

# include <math.h>
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>

# define MaxAtomNum   5000
# define MaxAtomNum2  5000
# define MaxAtomNum3  5000
# define N            10
# define Maxline      200
# define N1           15

FILE  *InpF;
FILE  *OutputFile;
char   OutputFName[Maxline], CoordFName[Maxline];
double xatom[MaxAtomNum], yatom[MaxAtomNum], zatom[MaxAtomNum];
char   atomname[MaxAtomNum][8];
char   atom_c_s[MaxAtomNum];
int    atomnumber, CoorFType;
double cell[N];
double TOdistmin, TOdistmax, OTOanglemin, OTOanglemax, TTdistmin, TTdistmax, TOTanglemin, TOTanglemax;
double TOdistave, OTOangleave, TOTangleave;
int    atomcalcnumber, atomcalc[MaxAtomNum];
int    countsingleT, countoxygen;
int    singleTatoms[5][MaxAtomNum2];
int    oxyg_neighb[3][MaxAtomNum3];
double atomcharge[MaxAtomNum];
int    siteTatoms[N1][N], siteOatoms[N1][N1]; // for T [1..12][1..5] for O [1..12][1..12]
int    numberTsites, numberOsites, ReferenceFType, numberGeatoms;
char   ReferenceFName[Maxline];
int    numberinTsite[N], numberinOsite[N1];   // for T [1..5] for O [1..12]
double xatomr[MaxAtomNum], yatomr[MaxAtomNum], zatomr[MaxAtomNum];
char   atomnamer[MaxAtomNum][8];
char   atom_c_sr[MaxAtomNum];
int    atom_input_order[MaxAtomNum], used_ref_atom[MaxAtomNum];


void   Initialization();
void   MainJob();
int    getatomcoor(char coordFName[Maxline], int CoorFType);
int    countwords(char line1[Maxline], int count1pos, int count2pos);
char  *grepnword(char *line1, int position_line, int nthword);
char  *grepnmword(char *line1, int position_ini, int position_final);
double atomsdist (double x1, double y1, double z1, double x2, double y2, double z2);
double atomsdistc(double x1, double y1, double z1, double x2, double y2, double z2);
double atomsangle(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33);
double atomsanglec(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33);
double absx(double x);
void   FindSingleT(char CoordFNamex[Maxline], int CoorFTypex);

double atomsangle_test(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33);

void main(argc, argv)
        int argc;
        char **argv;
{
int    crash_label, count1, count2, count3;
double dist, angle_degree;
double cryst_atom[4];  

    crash_label = 0.0;
    
    InpF = fopen(argv[1],"r");
   
    Initialization();
    
    FindSingleT(ReferenceFName, ReferenceFType);
    
    MainJob(); 
    
    printf(" thanks for using SiGe_TO_dist_angles_03 \n");
    
}     /******************************  end of main  **************************/

/******************************************************************************
function to initiate the job; reading general data and initialising
general variables.

******************************************************************************/
void    Initialization()
{   
    int  maxline, count1, count2;
    char line[Maxline];

    maxline = Maxline;     
    
    fgets(line, maxline, InpF);
    sscanf(line,"%s %d", CoordFName, &CoorFType);
    // name and type (gin, gout, xtl, car, cif ) of the input analyzing structure file 
    
    fgets(line, maxline, InpF);
    //sscanf(line,"%d %d %d", &numberGeatoms, &numberTsites, &numberOsites); 
    sscanf(line,"%d %d", &numberTsites, &numberOsites);
    // numberTsites, numberOsites
    
    fgets(line, maxline, InpF);
    sscanf(line,"%s %d", ReferenceFName, &ReferenceFType);
    // name and type (gin, gout, xtl, car, cif) of the input Reference pure silica structure file
    
    fgets(line, maxline, InpF);
    sscanf(line,"%lf %lf", &TOdistmin, &TOdistmax); 
    // min & max T-O distances 

    //fgets(line, maxline, InpF);
    //sscanf(line,"%s", OutputFName);
    // name of Output file 
      
    for (count1 = 1; count1 <= numberTsites; count1++)
    {
      fgets(line, maxline, InpF);
      sscanf(line,"%d", &numberinTsite[count1]);
      // multiplicity of site T [count1]
      
      for (count2 = 1; count2 <= numberinTsite[count1]; count2++)
      {
        fgets(line, maxline, InpF);
        sscanf(line,"%d", &siteTatoms[count2][count1]);
        // each T-atom [count2] at site T [count1]
        siteTatoms[count2][count1] = siteTatoms[count2][count1] - 1; // because later atoms arrays start from 0
      }
    }
    
    for (count1 = 1; count1 <= numberOsites; count1++)
    {
      fgets(line, maxline, InpF);
      sscanf(line,"%d", &numberinOsite[count1]);
      // multiplicity of site T [count1]
      
      for (count2 = 1; count2 <= numberinOsite[count1]; count2++)
      {
        fgets(line, maxline, InpF);
        sscanf(line,"%d", &siteOatoms[count2][count1]);
        // each T-atom [count2] at site T [count1]
        siteOatoms[count2][count1] = siteOatoms[count2][count1] - 1; // because later atoms arrays start from 0
      }
    }

    fclose(InpF);

    printf(" ms1 %s %d\n", CoordFName, CoorFType);
    //printf(" ms2 %d %d %d\n", numberGeatoms, numberTsites, numberOsites);
    printf(" ms2 %d %d\n", numberTsites, numberOsites);
    printf(" ms3 %s %d\n", ReferenceFName, ReferenceFType);
    printf(" ms4 %f %f\n", TOdistmin, TOdistmax);
    for (count1 = 1; count1 <= numberTsites; count1++)
    {
      printf(" ms5 %d T site and %d its elements below\n", count1, numberinTsite[count1]);
      
      for (count2 = 1; count2 <= numberinTsite[count1]; count2++)
      {
        printf("     %d  in Tsite_%d\n", siteTatoms[count2][count1], count1);
      }
    }
    for (count1 = 1; count1 <= numberOsites; count1++)
    {
      printf(" ms6 %d O site and %d its elements below\n", count1, numberinOsite[count1]);
      
      for (count2 = 1; count2 <= numberinOsite[count1]; count2++)
      {
        printf("     %d  in Osite_%d\n", siteOatoms[count2][count1], count1);
      }
    }
    
    printf("\n");
    
}    /*  end of Initialization   */

/******************************************************************************
function to do the Main Job

******************************************************************************/
void    MainJob()
{
    int    maxline, Maincount1, count1, count2, count3, count4, count5, count6, temp_char;
    char   line[Maxline], auxs1[Maxline], auxs2[Maxline], auxs3[Maxline], ginfile[Maxline];
    int    num_TOdist, num_OTOangle, num_TOTangle;
    double dist1, dist2, dist3, angle1, angle2;
    double value_runmin, value_runmax;
    long int lcount1;
    char   s1[Maxline], s2[Maxline], s3[Maxline], s4[Maxline];
    double tolerance_dist, disttmp;
    int    tmp_count1[N], tmp_count2[N1], tmp_count3, tmp_count4;
    int    go_on_1, go_on_2, go_on_3, Maincount2;
    double tolerance_delta;
    
    maxline = Maxline;
    
    //copy atom data to those of the reference
    strcpy(s1, ""); strcpy(s2, ""); strcpy(s3, ""); strcpy(s4, "");
    memset(s1, 0,  sizeof(s1)); memset(s2, 0,  sizeof(s2));
    memset(s3, 0,  sizeof(s3)); memset(s4, 0,  sizeof(s4));
    for (lcount1 = 0; lcount1 < MaxAtomNum; lcount1++)
    {
      //strcpy(atomname[lcount1], ""); memset(atomname[lcount1], 0,  sizeof(atomname[lcount1]));
      //strcpy(atomname[lcount1], ""); //memset(atomname, 0,  sizeof(atomname));
      for (count1 = 0; count1 < 8; count1++)
      {
        //atomname[lcount1][count1] = '\0';
        atomnamer[lcount1][count1] = '\0';
      }
      //strcpy(atom_c_s[count1], ""); memset(atom_c_s[count1], 0,  sizeof(atom_c_s[count1]));
    }
    for (count1 = 0; count1 < atomnumber; count1++)
    {
      xatomr[count1] = xatom[count1];
      yatomr[count1] = yatom[count1];
      zatomr[count1] = zatom[count1];
      strcpy(s1, atomname[count1]);
      sscanf(s1,"%s", atomnamer[count1]);
      strcpy(s3 , "c"); strcpy(s4 , "s");
      if (atom_c_s[count1] == 'c')
      {
        atom_c_sr[count1] = s3[0];
      }
      else
      {
        atom_c_sr[count1] = s4[0];
      }
      used_ref_atom[count1]    = 0;
      atom_input_order[count1] = -1;
    } 
    printf("Writing reference structure\n");
    for (count1 = 0; count1 <= atomnumber-1; count1++)
    {
      printf("%s %c %.5f %.5f %.5f\n", atomnamer[count1], atom_c_sr[count1], xatomr[count1], yatomr[count1], zatomr[count1]);  
    }
	printf("\n");    
    
    count1 = getatomcoor(CoordFName, CoorFType);
    
    //tolerance_dist = 0.45; tolerance_delta = 0.05; go_on_1 = 0;
    tolerance_dist = 1.30; tolerance_delta = 0.05; go_on_1 = 0;
    //while ( (tolerance_dist <= 1.30) && (go_on_1 == 0) ) 
    //{
      for (count1 = 0; count1 < atomnumber; count1++)
      {
        strcpy(s1, atomname[count1]);
        if ( ( (s1[0] =='s') || (s1[0] =='S') || (s1[0] =='A') || (s1[0] =='A') || (s1[0] =='G') || (s1[0] =='g') || (s1[0] =='o') || (s1[0] =='O') ) && (used_ref_atom[count1] == 0) )   
        {    
          for (count2 = 0; count2 < atomnumber; count2++)
          {
            strcpy(s2, atomname[count2]);
            if ( ( (s2[0] =='s') || (s2[0] =='S') || (s2[0] =='A') || (s2[0] =='A') || (s2[0] =='G') || (s2[0] =='g') || (s2[0] =='o') || (s2[0] =='O') ) && (used_ref_atom[count2] == 0) )
            {
              disttmp = atomsdist(xatomr[count2], yatomr[count2], zatomr[count2], xatom[count1], yatom[count1], zatom[count1]);
              if (disttmp <= tolerance_dist)
              {
                atom_input_order[count1] = count2; // gives the index of the atom in the reeference structure
                used_ref_atom[count1] = 1;
                printf("check reference 1  dist %f atom %3d atom_input_order %3d difference%3d\n", disttmp, count1, atom_input_order[count1], count1-atom_input_order[count1]);
                if ((count1-atom_input_order[count1]) != 0)
                {
                  printf("check reference 1a dist %f atom1 %3d %f %f %f atom2 %3d  %f %f %f\n", disttmp, count1, xatom[count1], yatom[count1], zatom[count1], count2, xatomr[count2], yatomr[count2]);
				}
                count2 = atomnumber + 1;
              }
            }
          }
        }
      }
      //tolerance_dist = tolerance_dist + tolerance_delta;
      count3 = 0;
      for (count4 = 0; count4 < atomnumber; count4++)
      {
        count3 = count3 + used_ref_atom[count4];
      }
      if (count3 == atomnumber)
      {
        go_on_1 = 1;
      }
    //}
    //if (go_on_1 == 1)
    //{
      //printf("GOOD atom order: reference and analyzing structures\n");
      for (count1 = 0; count1 < atomnumber; count1++)
      {
        //atom_input_order[count1] = count1;  // ojo esto es lo  que hayq quitar y corregir arriba adecuadamente
        printf("check reference 2 atom %3d atom_input_order %3d difference%3d\n", count1, atom_input_order[count1], count1-atom_input_order[count1]);
      }
    //}
    //else
    //{
    //  printf("BAD atom order: reference and analyzing structures\n");
    //  for (count1 = 0; count1 < atomnumber; count1++)
    //  {
    //    //atom_input_order[count1] = count1;  // ojo esto es lo  que hayq quitar y corregir arriba adecuadamente
    //    printf("%d %d %d\n", count1, atom_input_order[count1], count1-atom_input_order[count1]);
    //  }
    //}
    printf("\n");
    
    strncpy(OutputFName, CoordFName, strlen(CoordFName)-4);
    printf("CoordFName %s OutputFName %s\n", CoordFName, OutputFName);
    strcat(OutputFName, "_geom_data.txt");
    OutputFile = fopen(OutputFName,"w");

	fprintf(OutputFile, "Overall T (Si or Ge) Analysis\n");
    fprintf(OutputFile, "\n");

    // analyzing TO distances, overall
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
	  printf("atom %5d TO distances:", count1 + 1);
      for (count2 = 1; count2 <= 4; count2++)
      {
        count4 = singleTatoms[count2][Maincount1];
        count3 = atom_input_order[count4];
		if (count3 > -1)
        {
          dist1 = atomsdist(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
          if ( (dist1 < 1.3) || (dist1 > 2.1) )
          {
            printf(" OJO TO bad %f dist %.5f %.5f %.5f %.5f %.5f %.5f", dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
          }
          num_TOdist++;
          TOdistave = TOdistave + dist1;
          printf(" %f", dist1);
          if (dist1<value_runmin)
          {
            value_runmin = dist1;
          }
          if (dist1>value_runmax)
          {
            value_runmax = dist1;
          }
        }
      }
      printf("\n");  
    }  
    TOdistave = TOdistave / num_TOdist;
    fprintf(OutputFile, "TO_distances %f %f %f\n", value_runmin, TOdistave, value_runmax);
    // end of analyzing TO distances
    
    // analyzing OTO angles, overall
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
      if (Maincount1 == 27)
      {
        printf("Test error 0a atom %5d T %3d O1 %3d O2 %3d O3 %3d O4\n", singleTatoms[0][Maincount1], singleTatoms[1][Maincount1], singleTatoms[2][Maincount1], singleTatoms[3][Maincount1], singleTatoms[4][Maincount1]);
        count2 = singleTatoms[0][Maincount1]; count2 = atom_input_order[count2];
        printf("Test error 0b atom %5d T", count2);
        count2 = singleTatoms[1][Maincount1]; count2 = atom_input_order[count2];
        printf(" %3d O1", count2);
        count2 = singleTatoms[2][Maincount1]; count2 = atom_input_order[count2];
        printf(" %3d O2", count2);
        count2 = singleTatoms[3][Maincount1]; count2 = atom_input_order[count2];
        printf(" %3d O3", count2);
        count2 = singleTatoms[4][Maincount1]; count2 = atom_input_order[count2];
        printf(" %3d O4\n", count2);
      }
      printf("atom %5d at overall analysis OTO angles:", count1 + 1);
      for (count2 = 1; count2 < 4; count2++)
      {
      	tmp_count3 = singleTatoms[count2][Maincount1];
        count3     = atom_input_order[tmp_count3];
        if (count3 > -1)
        {
          for (count4 = count2 + 1; count4 <= 4; count4++)
          { 
            tmp_count4 = singleTatoms[count4][Maincount1];
            count5     = atom_input_order[tmp_count4];
            if (Maincount1 == 27)
            {
              printf("Test error 1a atom %5d T %3d O1 %3d O2\n", Maincount1+1, count3+1, count5+1);
              angle1 = atomsangle_test(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
            }
            if (count5 > -1)
            {
              angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);

              if ( (angle1 < 90.0) || (angle1 > 130.0) )
              {
                printf(" OJO OTO %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", angle1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
              }
              num_TOdist++;
              TOdistave = TOdistave + angle1;
              printf(" %f", angle1);
              if (angle1<value_runmin)
              {
                value_runmin = angle1;
              }
              if (angle1>value_runmax)
              {
                value_runmax = angle1;
              }
            }
            else
            {
              printf("Test error 2b atom %5d bad labelled oxygen neighbour %3d giving %3d at first oxygen %3d\n", Maincount1, count4, count5, count2);
            }
          }
        }
        else
        {
          printf("Test error 2a atom %5d bad labelled oxygen neighbour %3d giving %3d\n", Maincount1, count2, count3);
        }
      }
      printf("\n");  
    }  
    TOdistave = TOdistave / num_TOdist;
    fprintf(OutputFile, "OTO_angles %f %f %f\n", value_runmin, TOdistave, value_runmax);
    // end of analyzing OTO angles
    
    // analyzing TOT angles, overall
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      if ( (count5 > -1) && (count3 > -1) )
      {
        angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
   
        if ( (angle1 < 110.0) || (angle1 >= 178.0) )
        {
          printf(" atom %5d TOT angle: OJO TOT %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", Maincount1, dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        num_TOdist++;
        TOdistave = TOdistave + angle1;
        printf("atom %5d TOT angle: %f\n", Maincount1, angle1);
        if (angle1<value_runmin)
        {
          value_runmin = angle1;
        }
        if (angle1>value_runmax)
        {
          value_runmax = angle1;
        }
      }
    }  
    TOdistave = TOdistave / num_TOdist;
    fprintf(OutputFile, "TOT_angles %f %f %f\n", value_runmin, TOdistave, value_runmax);
    if (value_runmax>=175.0)
    {
      printf("OJO Linear Angle See Careful %f\n", value_runmax);
    }
    // end of analyzing TOT angles
    
    // analyzing TT distances, overall
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      if ( (count5 > -1) && (count3 > -1) )
      {
        dist1 = atomsdist(xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        num_TOdist++;
        TOdistave = TOdistave + dist1;
        if ( (dist1 < 2.5) || (dist1 > 4) )
        {
          printf("atom %5d TT distance: %f OJO TT bad %.5f %.5f %.5f %.5f %.5f %.5f\n", Maincount1, dist1, xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        else
        {
          printf("atom %5d TT distance: %f\n", Maincount1, dist1);
        }
        if (dist1<value_runmin)
        {
          value_runmin = dist1;
        }
        if (dist1>value_runmax)
        {
          value_runmax = dist1;
        }
      }
    }  
    TOdistave = TOdistave / num_TOdist;
    fprintf(OutputFile, "TT_distances %f %f %f\n", value_runmin, TOdistave, value_runmax);
    // end of analyzing TT angles
    
    // here is the segment for overall focalizing Ge containing analysis
    fprintf(OutputFile, "\n");
	fprintf(OutputFile, "Ge Analysis\n");
    fprintf(OutputFile, "\n");
    // analyzing TO distances, Focalized on Ge
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
	  printf("atom %5d TO distances:", count1 + 1);
      strcpy(s1, atomname[count1]);
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        for (count2 = 1; count2 <= 4; count2++)
        {
          count4 = singleTatoms[count2][Maincount1];
          count3 = atom_input_order[count4];
		  if (count3 > -1)
          {
            dist1 = atomsdist(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
            if ( (dist1 < 1.3) || (dist1 > 2.1) )
            {
              printf(" OJO TO bad %f dist %.5f %.5f %.5f %.5f %.5f %.5f", dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
            }
            num_TOdist++;
            TOdistave = TOdistave + dist1;
            printf(" %f", dist1);
            if (dist1<value_runmin)
            {
              value_runmin = dist1;
            }
            if (dist1>value_runmax)
            {
              value_runmax = dist1;
            }
          }
        }
      }
      printf("\n");  
    }
	if (num_TOdist == 0)
	{
      fprintf(OutputFile, "GeO_distances 0.0 0.0 0.0\n");
	}
	else
	{
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "GeO_distances %f %f %f\n", value_runmin, TOdistave, value_runmax);
    }
    // end of analyzing TO distances, Focalized on Ge
    
    // analyzing OTO angles, Focalized on Ge
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
      printf("atom %5d OTO angles:", count1 + 1);
      strcpy(s1, atomname[count1]);
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        for (count2 = 1; count2 < 4; count2++)
        {
          tmp_count3 = singleTatoms[count2][Maincount1];
          count3     = atom_input_order[tmp_count3];
          if (count3 > -1)
          {
            for (count4 = count2 + 1; count4 <= 4; count4++)
            { 
              tmp_count4 = singleTatoms[count4][Maincount1];
              count5     = atom_input_order[tmp_count4];
              if (count5 > -1)
              {
                angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
  
                if ( (angle1 < 90.0) || (angle1 > 130.0) )
                {
                  printf(" OJO OTO %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
                }
                num_TOdist++;
                TOdistave = TOdistave + angle1;
                printf(" %f", angle1);
                if (angle1<value_runmin)
                {
                  value_runmin = angle1;
                }
                if (angle1>value_runmax)
                {
                  value_runmax = angle1;
                }
              }
            }
          }
        }
        printf("\n");  
      }
    } 
    if (num_TOdist == 0)
	{
      fprintf(OutputFile, "OGeO_angles 0.0 0.0 0.0\n");
	}
	else
    {
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "OGeO_angles %f %f %f\n", value_runmin, TOdistave, value_runmax);
    }
    // end of analyzing OTO angles, Focalized on Ge
    
    // analyzing TOT angles, Focalized on Ge
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      strcpy(s1, atomname[count3]); tmp_count3 = 0;
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        tmp_count3 = 1;
      }  
      count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      strcpy(s2, atomname[count5]); tmp_count4 = 0;
      if ( (s2[0] =='g') || (s2[0] =='G') )
      {
        tmp_count4 = 1;
      }
      if ( (count5 > -1) && (count3 > -1) && ( (tmp_count3 == 1) || (tmp_count4 == 1) ) )
      {
        angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
   
        if ( (angle1 < 110.0) || (angle1 >= 178.0) )
        {
          printf(" atom %5d TOT angle: OJO TOT %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", Maincount1, dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        num_TOdist++;
        TOdistave = TOdistave + angle1;
        printf("atom %5d TOT angle: %f\n", Maincount1, angle1);
        if (angle1<value_runmin)
        {
          value_runmin = angle1;
        }
        if (angle1>value_runmax)
        {
          value_runmax = angle1;
        }
      }
    }  
    if (num_TOdist == 0)
	{
      fprintf(OutputFile, "TOT_angles 0.0 0.0 0.0 (T1=Ge, T2=Ge o Si)\n");
	}
	else
	{
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "TOT_angles %f %f %f (T1=Ge, T2=Ge o Si)\n", value_runmin, TOdistave, value_runmax);
    }
    if (value_runmax>=175.0)
    {
      printf("OJO Linear Angle See Careful %f\n", value_runmax);
    }
    // end of analyzing TOT angles, Focalized on Ge
    
    // analyzing TT distances, Focalized on Ge
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      strcpy(s1, atomname[count3]); tmp_count3 = 0;
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        tmp_count3 = 1;
      } 
      count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      strcpy(s2, atomname[count5]); tmp_count4 = 0;
      if ( (s2[0] =='g') || (s2[0] =='G') )
      {
        tmp_count4 = 1;
      }
      if ( (count5 > -1) && (count3 > -1) && ( (tmp_count3 == 1) || (tmp_count4 == 1) ) )
      {
        dist1 = atomsdist(xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        num_TOdist++;
        TOdistave = TOdistave + dist1;
        if ( (dist1 < 2.5) || (dist1 > 4) )
        {
          printf("atom %5d TT distance: %f OJO TT bad %.5f %.5f %.5f %.5f %.5f %.5f\n", Maincount1, dist1, xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        else
        {
          printf("atom %5d TT distance: %f\n", Maincount1, dist1);
        }
        if (dist1<value_runmin)
        {
          value_runmin = dist1;
        }
        if (dist1>value_runmax)
        {
          value_runmax = dist1;
        }
      }
    }  
    if (num_TOdist == 0)
	{
      fprintf(OutputFile, "TT_distances 0.0 0.0 0.0 (T1=Ge, T2=Ge o Si)\n");
	}
	else
	{
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "TT_distances %f %f %f (T1=Ge, T2=Ge o Si)\n", value_runmin, TOdistave, value_runmax);
    }
    // end of analyzing TT angles, Focalized on Ge
    
    // end of the segment for overall focalizing Ge containing analysis
    
    // segment of overall analysis by T specific sites
	fprintf(OutputFile, "\n");
	fprintf(OutputFile, "Overall T specific T (Si or Ge) Analysis\n");
    fprintf(OutputFile, "\n");

    // analyzing TO distances, T specific
  for (Maincount2 = 1; Maincount2 <= numberTsites; Maincount2++)
  {
  	//fprintf(OutputFile, "T %d specific\n", Maincount2);
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
	  printf("atom %5d %d of site TO distances:", count1 + 1, Maincount2);
	  go_on_1 = 0;
	  for (count2 = 1; count2 <= numberinTsite[Maincount2]; count2++)
      {
        if (count1 == siteTatoms[count2][Maincount2])
        {
          go_on_1 = 1; 
          count2 <= numberinTsite[Maincount2] + 10;  
        }
      }
	  if (go_on_1 == 1)
      {  
        for (count2 = 1; count2 <= 4; count2++)
        {
          count4 = singleTatoms[count2][Maincount1];
          count3 = atom_input_order[count4];
		  if (count3 > -1)
          {
            dist1 = atomsdist(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
            if ( (dist1 < 1.3) || (dist1 > 2.1) )
            {
              printf(" OJO TO bad %f dist %.5f %.5f %.5f %.5f %.5f %.5f", dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
            }
            num_TOdist++;
            TOdistave = TOdistave + dist1;
            printf(" %f", dist1);
            if (dist1<value_runmin)
            {
              value_runmin = dist1;
            }
            if (dist1>value_runmax)
            {
              value_runmax = dist1;
            }
          }
        }
        printf("\n");  
      }
    }
	if (num_TOdist == 0)
	{
      fprintf(OutputFile, "TO_distances T %d specific 0.0 0.0 0.0\n", Maincount2);
	}
	else
	{  
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "TO_distances T %d specific %f %f %f\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }   
    // end of analyzing TO distances
  }
  
  for (Maincount2 = 1; Maincount2 <= numberTsites; Maincount2++)
  {
    // analyzing OTO angles, T specific
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
      printf("atom %5d OTO angles:", count1 + 1);
      go_on_1 = 0;
	  for (count2 = 1; count2 <= numberinTsite[Maincount2]; count2++)
      {
        if (count1 == siteTatoms[count2][Maincount2])
        {
          go_on_1 = 1; 
          count2 <= numberinTsite[Maincount2] + 10;  
        }
      }
	  if (go_on_1 == 1)
      {
      	
        for (count2 = 1; count2 < 4; count2++)
        {
          tmp_count3 = singleTatoms[count2][Maincount1];
          count3     = atom_input_order[tmp_count3];
          if (count3 > -1)
          {
            for (count4 = count2 + 1; count4 <= 4; count4++)
            { 
              tmp_count4 = singleTatoms[count4][Maincount1];
              count5     = atom_input_order[tmp_count4];
              if (count5 > -1)
              {
                angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
  
                if ( (angle1 < 90.0) || (angle1 > 130.0) )
                {
                  printf(" OJO OTO %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
                }
                num_TOdist++;
                TOdistave = TOdistave + angle1;
                printf(" %f", angle1);
                if (angle1<value_runmin)
                {
                  value_runmin = angle1;
                }
                if (angle1>value_runmax)
                {
                  value_runmax = angle1;
                }
              }
            }
          }
        }
        printf("\n");  
      }
    }
	if (num_TOdist == 0)
	{
      fprintf(OutputFile, "OTO_angles T %d specific 0.0 0.0 0.0\n", Maincount2);
	}
	else
	{  
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "OTO_angles T %d specific %f %f %f\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }
    // end of analyzing OTO angles
  }
  
  for (Maincount2 = 1; Maincount2 <= numberOsites; Maincount2++)
  {
    // analyzing TOT angles, T specific
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      go_on_1 = 0;
      for (count2 = 1; count2 <= numberinOsite[Maincount2]; count2++)
      {
        if (count1 == siteOatoms[count2][Maincount2])
        {
          go_on_1 = 1; 
          count2 <= numberinOsite[Maincount2] + 10;  
        }
      }
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      if ( (count5 > -1) && (count3 > -1) && (go_on_1 == 1) )
      {
        angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
   
        if ( (angle1 < 110.0) || (angle1 >= 178.0) )
        {
          printf(" atom %5d TOT angle: OJO TOT %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", Maincount1, dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        num_TOdist++;
        TOdistave = TOdistave + angle1;
        printf("atom %5d TOT angle: %f\n", Maincount1, angle1);
        if (angle1<value_runmin)
        {
          value_runmin = angle1;
        }
        if (angle1>value_runmax)
        {
          value_runmax = angle1;
        }
      }
    }  
    if (num_TOdist == 0)
	{
      fprintf(OutputFile, "TOT_angles O %d specific 0.0 0.0 0.0\n", Maincount2);
	}
	else
	{
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "TOT_angles O %d specific %f %f %f\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }
    if (value_runmax>=175.0)
    {
      printf("OJO Linear Angle See Careful %f\n", value_runmax);
    }
    // end of analyzing TOT angles
  }
  
  for (Maincount2 = 1; Maincount2 <= numberOsites; Maincount2++)
  {
    // analyzing TT distances, T specific
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      go_on_1 = 0;
      for (count2 = 1; count2 <= numberinOsite[Maincount2]; count2++)
      {
        if (count1 == siteOatoms[count2][Maincount2])
        {
          go_on_1 = 1; 
          count2 <= numberinOsite[Maincount2] + 10;  
        }
      }
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      if ( (count5 > -1) && (count3 > -1) && (go_on_1 == 1) )
      {
        dist1 = atomsdist(xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        num_TOdist++;
        TOdistave = TOdistave + dist1;
        if ( (dist1 < 2.5) || (dist1 > 4) )
        {
          printf("atom %5d TT distance: %f OJO TT bad %.5f %.5f %.5f %.5f %.5f %.5f\n", Maincount1, dist1, xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        else
        {
          printf("atom %5d TT distance: %f\n", Maincount1, dist1);
        }
        if (dist1<value_runmin)
        {
          value_runmin = dist1;
        }
        if (dist1>value_runmax)
        {
          value_runmax = dist1;
        }
      }
    }  
    if (num_TOdist == 0)
	{
      fprintf(OutputFile, "TT_distances O %d specific 0.0 0.0 0.0\n", Maincount2);
	}
	else
	{
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "TT_distances O %d specific %f %f %f\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }
    // end of analyzing TT angles
  }	
	// end of overall analysis by T specific sites 
	
	    // segment of Ge-particular analysis by T specific sites
	fprintf(OutputFile, "\n");
	fprintf(OutputFile, "Ge-particular T specific T (Si or Ge) Analysis\n");
    fprintf(OutputFile, "\n");

    // analyzing TO distances, T specific
  for (Maincount2 = 1; Maincount2 <= numberTsites; Maincount2++)
  {
  	//fprintf(OutputFile, "T %d specific\n", Maincount2);
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
	  printf("atom %5d %d of site TO distances:", count1 + 1, Maincount2);
	  go_on_1 = 0;
	  strcpy(s1, atomname[count1]);
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
	    for (count2 = 1; count2 <= numberinTsite[Maincount2]; count2++)
        {
          if (count1 == siteTatoms[count2][Maincount2])
          {
            go_on_1 = 1; 
            count2 <= numberinTsite[Maincount2] + 10;  
          }
        }
      }
	  if (go_on_1 == 1)
      {  
        for (count2 = 1; count2 <= 4; count2++)
        {
          count4 = singleTatoms[count2][Maincount1];
          count3 = atom_input_order[count4];
		  if (count3 > -1)
          {
            dist1 = atomsdist(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
            if ( (dist1 < 1.3) || (dist1 > 2.1) )
            {
              printf(" OJO TO bad %f dist %.5f %.5f %.5f %.5f %.5f %.5f", dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
            }
            num_TOdist++;
            TOdistave = TOdistave + dist1;
            printf(" %f", dist1);
            if (dist1<value_runmin)
            {
              value_runmin = dist1;
            }
            if (dist1>value_runmax)
            {
              value_runmax = dist1;
            }
          }
        }
        printf("\n");  
      }
    }
	if (num_TOdist == 0)
	{
      fprintf(OutputFile, "GeO_distances Ge %d specific 0.0 0.0 0.0\n", Maincount2);
	}
	else
	{  
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "GeO_distances Ge %d specific %f %f %f\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }   
    // end of analyzing TO distances
  }
  
  for (Maincount2 = 1; Maincount2 <= numberTsites; Maincount2++)
  {
    // analyzing OTO angles, T specific
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count2 = singleTatoms[0][Maincount1];
      count1 = atom_input_order[count2];
      printf("atom %5d OTO angles:", count1 + 1);
      go_on_1 = 0;
      strcpy(s1, atomname[count1]); 
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
	    for (count2 = 1; count2 <= numberinTsite[Maincount2]; count2++)
        {
          if (count1 == siteTatoms[count2][Maincount2])
          {
            go_on_1 = 1; 
            count2 <= numberinTsite[Maincount2] + 10;  
          }
        }
      }   
	  if (go_on_1 == 1)
      {
      	
        for (count2 = 1; count2 < 4; count2++)
        {
          tmp_count3 = singleTatoms[count2][Maincount1];
          count3     = atom_input_order[tmp_count3];
          if (count3 > -1)
          {
            for (count4 = count2 + 1; count4 <= 4; count4++)
            { 
              tmp_count4 = singleTatoms[count4][Maincount1];
              count5     = atom_input_order[tmp_count4];
              if (count5 > -1)
              {
                angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
  
                if ( (angle1 < 90.0) || (angle1 > 130.0) )
                {
                  printf(" OJO OTO %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
                }
                num_TOdist++;
                TOdistave = TOdistave + angle1;
                printf(" %f", angle1);
                if (angle1<value_runmin)
                {
                  value_runmin = angle1;
                }
                if (angle1>value_runmax)
                {
                  value_runmax = angle1;
                }
              }
            }
          }
        }
        printf("\n");  
      }
    }
	if (num_TOdist == 0)
	{
      fprintf(OutputFile, "OGeO_angles Ge %d specific 0.0 0.0 0.0\n", Maincount2);
	}
	else
	{  
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "OGeO_angles Ge %d specific %f %f %f\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }
    // end of analyzing OTO angles
  }
  
  for (Maincount2 = 1; Maincount2 <= numberOsites; Maincount2++)
  {
    // analyzing TOT angles, T specific
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      go_on_1 = 0;
      for (count2 = 1; count2 <= numberinOsite[Maincount2]; count2++)
      {
        if (count1 == siteOatoms[count2][Maincount2])
        {
          go_on_1 = 1; 
          count2 <= numberinOsite[Maincount2] + 10;  
        }
      }
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      strcpy(s1, atomname[count3]); tmp_count3 = 0;
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        tmp_count3 = 1;
      }
	  count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      strcpy(s1, atomname[count5]); tmp_count4 = 0;
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        tmp_count4 = 1;
      }
      if ( (count5 > -1) && (count3 > -1) && (go_on_1 == 1) && ( (tmp_count3 == 1) || (tmp_count4 == 1) ) )
      {
        angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
   
        if ( (angle1 < 110.0) || (angle1 >= 178.0) )
        {
          printf(" atom %5d TOT angle: OJO TOT %f bad %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", Maincount1, dist1, xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        num_TOdist++;
        TOdistave = TOdistave + angle1;
        printf("atom %5d TOT angle: %f\n", Maincount1, angle1);
        if (angle1<value_runmin)
        {
          value_runmin = angle1;
        }
        if (angle1>value_runmax)
        {
          value_runmax = angle1;
        }
      }
    }  
    if (num_TOdist == 0)
	{
      fprintf(OutputFile, "TOT_angles O %d specific 0.0 0.0 0.0 (T1=Ge, T2=Ge o Si)\n", Maincount2);
	}
	else
	{
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "TOT_angles O %d specific %f %f %f (T1=Ge, T2=Ge o Si)\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }
    if (value_runmax>=175.0)
    {
      printf("OJO Linear Angle See Careful %f\n", value_runmax);
    }
    // end of analyzing TOT angles
  }
  
  for (Maincount2 = 1; Maincount2 <= numberOsites; Maincount2++)
  {
    // analyzing TT distances, T specific
    num_TOdist = 0; TOdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countoxygen; Maincount1++)
    {
      count2 = oxyg_neighb[0][Maincount1];
      count1 = atom_input_order[count2];
      go_on_1 = 0;
      for (count2 = 1; count2 <= numberinOsite[Maincount2]; count2++)
      {
        if (count1 == siteOatoms[count2][Maincount2])
        {
          go_on_1 = 1; 
          count2 <= numberinOsite[Maincount2] + 10;  
        }
      }
      count4 = oxyg_neighb[1][Maincount1];
      count3 = atom_input_order[count4];
      strcpy(s1, atomname[count3]); tmp_count3 = 0;
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        tmp_count3 = 1;
      }
      count6 = oxyg_neighb[2][Maincount1];
      count5 = atom_input_order[count6];
      strcpy(s1, atomname[count5]); tmp_count4 = 0;
      if ( (s1[0] =='g') || (s1[0] =='G') )
      {
        tmp_count4 = 1;
      }
      if ( (count5 > -1) && (count3 > -1) && (go_on_1 == 1) && ( (tmp_count3 == 1) || (tmp_count4 == 1) ) )
      {
        dist1 = atomsdist(xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        num_TOdist++;
        TOdistave = TOdistave + dist1;
        if ( (dist1 < 2.5) || (dist1 > 4) )
        {
          printf("atom %5d TT distance: %f OJO TT bad %.5f %.5f %.5f %.5f %.5f %.5f\n", Maincount1, dist1, xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        }
        else
        {
          printf("atom %5d TT distance: %f\n", Maincount1, dist1);
        }
        if (dist1<value_runmin)
        {
          value_runmin = dist1;
        }
        if (dist1>value_runmax)
        {
          value_runmax = dist1;
        }
      }
    }  
    if (num_TOdist == 0)
	{
      fprintf(OutputFile, "TT_distances O %d specific 0.0 0.0 0.0 (T1=Ge, T2=Ge o Si)\n", Maincount2);
	}
	else
	{
      TOdistave = TOdistave / num_TOdist;
      fprintf(OutputFile, "TT_distances O %d specific %f %f %f (T1=Ge, T2=Ge o Si)\n", Maincount2, value_runmin, TOdistave, value_runmax);
    }
    // end of analyzing TT angles
  }	
	// end of Ge-particular analysis by T specific sites 

    
    fclose(OutputFile);

}  /************* end of MainJob *************/

/****************************************************************************
     this function extractes the atoms names, core-shel specification
     and coordinates

     types of files allowed
     1:   gulp input and/or restar   *.gin, *.res
     2:   gulp output                *.gout or *.gt
     3:   msi xtl                    *.xtl
     4:   gulp car file              *.car  NOTE, is GULP CAR FILE
     5:   cif                        *.cif

     the results are storaged in external variables:
     double xatom[MaxAtomNum], yatom[MaxAtomNum], zatom[MaxAtomNum];
     char   atomname[MaxAtomNum][4];
     char   atom_c_s[MaxAtomNum];
     int    atomnumber;
     double cell[N];
*****************************************************************************/
int    getatomcoor(char coordFName[Maxline], int CoorFType)
{
        double r1, r2, r3, r4, tempval1;
        int    c1, c2, c3, c4, tmp_dummy, cont1, maxline, intspace, contspace;
        char   s1[Maxline], s2[Maxline], s3[Maxline], s4[Maxline];
        char   s5[Maxline], s6[Maxline], s7[Maxline], s8[Maxline];
        char   line[Maxline];
        int    readinline;
        int    getatomcoor;
        FILE  *InProcesFile;
        int    count1, count2, count3;
        char  *chtmp1;
        char   s1a[Maxline];
        long int lcount1;

        // initializing string variables
        strcpy(s1, ""); strcpy(s2, ""); strcpy(s3, ""); strcpy(s4, "");
        strcpy(s5, ""); strcpy(s6, ""); strcpy(s7, ""); strcpy(s8, "");
        memset(s1, 0,  sizeof(s1)); memset(s1, 0,  sizeof(s2));
        memset(s1, 0,  sizeof(s3)); memset(s1, 0,  sizeof(s4));
        memset(s1, 0,  sizeof(s5)); memset(s1, 0,  sizeof(s6));
        memset(s1, 0,  sizeof(s7)); memset(s1, 0,  sizeof(s8));
        strcpy(line, ""); memset(s1, 0,  sizeof(line));
        strcpy(s1a, ""); memset(s1a, 0,  sizeof(s1a));
        for (lcount1 = 0; lcount1 < MaxAtomNum; lcount1++)
        {
          //strcpy(atomname[lcount1], ""); memset(atomname[lcount1], 0,  sizeof(atomname[lcount1]));
          //strcpy(atomname[lcount1], ""); //memset(atomname, 0,  sizeof(atomname));
          for (count1 = 0; count1 < 8; count1++)
          {
            //atomname[lcount1][count1] = '\0';
            atomname[lcount1][count1] = '\0';
          }
//          strcpy(atom_c_s[count1], ""); memset(atom_c_s[count1], 0,  sizeof(atom_c_s[count1]));
        }


        printf("already in getatomcoor 0, reading file %s of type %d \n", coordFName, CoorFType);
        
        InProcesFile = fopen(coordFName,"r");
        getatomcoor = 1; 
        maxline    = Maxline;
        atomnumber = 0;
        readinline = 1;
        
        printf("already in getatomcoor 1, reading file %s of type %d \n", coordFName, CoorFType);
        
        switch (CoorFType)
        {
        case 1:
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "cell", 4)!= 0)
           {
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           fgets(line, maxline, InProcesFile);
           /* blanck line skipper */
           if (line == NULL)
             {
             getatomcoor = 0;
             break;
             }
           while (blanckline(line) == 0)
              {
                fgets(line, maxline, InProcesFile);
              }
           /*  end of blanck line skipper  */
           sscanf(line,"%lf %lf %lf %lf %lf %lf", &cell[1], &cell[2], &cell[3], &cell[4], &cell[5], &cell[6]);
           strcpy(s1, "fractional");
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
           while (readinline == 1)
             {
               fgets(line, maxline, InProcesFile);
               if (line == NULL)
                {
                  getatomcoor = 0;
                  break;
                }
               while (blanckline(line) == 0)
                {
                  fgets(line, maxline, InProcesFile);
                }

/*               printf(" aviso0 %s\n", line);      */
               if (line != NULL)
                 {
                   if ( countwords(line, 0,strlen(line)-1) <= 3)
                     {
                       readinline = 0;
                     }
                     else
                     {
/*                       printf(" aviso1 %d\n", readinline);     */
                       tmp_dummy = 1;
                       //printf(" palabra1a %s\n", s2);
                       //strcpy(s2 , grepnword(line, 0, 1));
                       sscanf(line,"%s", s2);
                       //printf(" palabra1b %s\n", s2);           
                       strcpy(s3 , grepnword(line, 0, 2));
/*                       printf(" palabra2 %s\n", s3);           */
                       strcpy(s4 , grepnword(line, 0, 3));
/*                       printf(" palabra3 %s\n", s4);           */
                       strcpy(s5 , grepnword(line, 0, 4));
/*                       printf(" palabra4 %s\n", s5);           */
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            if ( countwords(line, 0,strlen(line)-1) <= 4)
                              readinline = 0;
                            else
                              strcpy(s6 , grepnword(line, 0, 5));
/*                            printf(" palabra5 %s\n", s6);      */
                            if (readinline != 0)
                               {
/*                                 printf(" aviso s6a %d\n", readinline);   */
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     {
                                       tmp_dummy = 1;
                                     }
                                   else
                                     {
                                       readinline = 0;
                                     }
                               }
                               else
                                 {
                                   readinline = 0;
/*                                   printf(" aviso s6b\n");      */
                                 }
                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {  
                                 readinline = 0;
                               }
                               
                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {
                                 readinline = 0;
                               }
                                                          
/*                            printf(" aviso readinline = %d\n", readinline);       */
                            if (readinline != 0)
                              {
                                atomnumber++;
                                sscanf(s4, "%lf", &xatom[atomnumber-1]);
/*                                printf("%f\t", xatom[atomnumber-1]);   */
                                sscanf(s5, "%lf", &yatom[atomnumber-1]);
/*                                printf(" %f\t", yatom[atomnumber-1]);  */
                                sscanf(s6, "%lf", &zatom[atomnumber-1]);
/*                                printf("%f\n", zatom[atomnumber-1]);   */
                                
//                                tmp_dummy = 0;
//                                for (tmp_dummy = 0; tmp_dummy <= strlen(s2) -1; tmp_dummy++)
//                                   atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy];                                
                                sscanf(s2,"%s", atomname[atomnumber-1]); 
                                /*
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                */
                                tmp_dummy = 1;
/*                                printf(" nombre atomo ?? %d %s\n", atomnumber-1, atomname[atomnumber-1]);  */
                               
                                atom_c_s[atomnumber-1] = s3[0];
                              }
                          }
                          else
                          {
                            printf(" llega a aqui ???");            
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                sscanf(s2,"%s", atomname[atomnumber-1]);
//                                tmp_dummy = 1;
//                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
//                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
/*                          printf(" retorno de la columna 27\n");  */
                          }
/*                     printf(" retorno de la columna 22\n");       */
                     }
/*                 printf(" retorno de la columna 18\n");           */
                 }
                 else
                   readinline = 0;  
/*             printf(" retorno de la columna 14\n");               */
             }
             
             printf("already in getatomcoor 2\n");

           break;

        case 2:
           
           strcpy(s1, "Final fractional coordinates of atoms");
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
           for (c1 = 1; c1 <= 6; c1++)
               {
               fgets(line, maxline, InProcesFile);
               }
           while (readinline == 1)
             {
               sscanf(line,"%s", s1);
               while ( (line != NULL) && (isspace(s1[0]) == 0) )
                  {
                  fgets(line, maxline, InProcesFile);
                  sscanf(line,"%s", s1);
                  }

               if (line != NULL)
                 {
                   if ( (strncmp(line, "--------------------", 20) != 0) || (strncmp(line, "Final cell", 10)!= 0) )
                     {
                       readinline = 0;
                     }
                     else
                     {
                       tmp_dummy = 1;
                       strcpy(s2 , grepnword(line, 0, 2));
                       strcpy(s3 , grepnword(line, 0, 3));
                       strcpy(s4 , grepnword(line, 0, 4));
                       strcpy(s5 , grepnword(line, 0, 5));
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            strcpy(s6 , grepnword(line, 0, 6));
                            if (strstr(s6, "") == NULL)
                               {
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     tmp_dummy = 1;
                                   else
                                     readinline = 0;
                               }
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf(s3[0], "%c", atom_c_s[atomnumber-1]);
                                sprintf(s4, "%.5f",xatom[atomnumber-1]);
                                sprintf(s5, "%.5f",yatom[atomnumber-1]);
                                sprintf(s6, "%.5f",zatom[atomnumber-1]);
                              }
                          }
                          else
                          {
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
                          }
                     }
                 }
                 else
                   readinline = 0;  
             }
           
           strcpy(s1, "Final cell parameters and derivatives");
           tmp_dummy = 1000;
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
                                      
           if (tmp_dummy != 0)
           {
             for (c1 = 1; c1 <= 6; c1++)
                 cell[c1] = -1.0;
                 c3 = 1; 
           }                                                                                                                                       
           else
           {
             for (c1 = 1; c1 <= 3; c1++)
               {
               fgets(line, maxline, InProcesFile);
               }
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[1], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s", s1, &cell[2], s2);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[3], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[4], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[5], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[6], s2, s3, &tempval1, s4);
           }
           break;

        case 3:
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "CELL", 4)!= 0)
           {   
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           fgets(line, maxline, InProcesFile);
           /* blanck line skipper */
           if (line == NULL)
             {
             getatomcoor = 0;
             break;
             }
           while (blanckline(line) == 0)
              {
                fgets(line, maxline, InProcesFile);
              }
           /*  end of blanck line skipper  */
           sscanf(line,"%lf %lf %lf %lf %lf %lf", &cell[1], &cell[2], &cell[3], &cell[4], &cell[5], &cell[6]);
           strcpy(s1, "NAME");
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
/*           fgets(line, maxline, InProcesFile);   */
           while (readinline == 1)
             {
               fgets(line, maxline, InProcesFile);
               if (line == NULL)
                {
                  getatomcoor = 0;
                  break;
                }
               while (blanckline(line) == 0)
                {
                  fgets(line, maxline, InProcesFile);
                }

/*               printf(" aviso0 %s\n", line);     */
               if (line != NULL)
                 {
                   if ( countwords(line, 0,strlen(line)-1) <= 3)
                     {
                       readinline = 0;
                     }
                     else
                     {
/*                       printf(" aviso1 %d\n", readinline);     
                       tmp_dummy = 1;
                       strcpy(s2 , grepnword(line, 0, 1));
                      printf(" palabra1 %s\n", s2);           
                       strcpy(s3 , "c");
                       printf(" palabra2 %s\n", s3);           
                       strcpy(s4 , grepnword(line, 0, 2));
                       printf(" palabra3 %s\n", s4);           
                       strcpy(s5 , grepnword(line, 0, 3));
*/
                       sscanf(line,"%s %s %s %s", s2, s4, s5, s6);
                       strcpy(s3 , "c");
/*                       printf(" read %s %s %s %s\n", s2, s3, s4, s5);          */
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            if ( countwords(line, 0,strlen(line)-1) <= 3)
                              readinline = 0;
                            else
                               tmp_dummy = 1;
/*                              strcpy(s6 , grepnword(line, 0, 4)); */
/*                            printf(" palabra5 %s\n", s6);      */
                            if (readinline != 0)
                               {
/*                                 printf(" aviso s6a %d\n", readinline);   */
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     {
                                       tmp_dummy = 1;
                                     }
                                   else
                                     {
                                       readinline = 0;
                                     }
                               }
                               else
                                 {
                                   readinline = 0;
/*                                   printf(" aviso s6b\n");      */
                                 }
                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {  
                                 readinline = 0;
                               }
                               
                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {
                                 readinline = 0;
                               }
                                                          
/*                            printf(" aviso readinline = %d\n", readinline);       */
                            if (readinline != 0)
                              {
                                atomnumber++;
                                sscanf(s4, "%lf", &xatom[atomnumber-1]);
/*                                printf("%f\t", xatom[atomnumber-1]);   */
                                sscanf(s5, "%lf", &yatom[atomnumber-1]);
/*                                printf(" %f\t", yatom[atomnumber-1]);  */
                                sscanf(s6, "%lf", &zatom[atomnumber-1]);
/*                                printf("%f\n", zatom[atomnumber-1]);   */
                                
                                tmp_dummy = 0;
                                for (tmp_dummy = 0; tmp_dummy <= strlen(s2) -1; tmp_dummy++)
                                   atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy];                                
                                /*
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                */
                                tmp_dummy = 1;
/*                                printf(" nombre atomo ?? %d %s\n", atomnumber-1, atomname[atomnumber-1]);  */
                               
                                atom_c_s[atomnumber-1] = s3[0];
/*                                printf(" %c\n", atom_c_s[atomnumber-1]);                               
                                                              
                                printf("ojo %d\t %s\t %c\t", atomnumber-1,  atomname[atomnumber-1], atom_c_s[atomnumber-1]);
                                printf("%f\t %f\t %f; %s\n", xatom[atomnumber-1], yatom[atomnumber-1], zatom[atomnumber-1], line); 
*/
                              }
                          }
                          else
                          {
                            printf(" llega a aqui ???");            
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
/*                          printf(" retorno de la columna 27\n");  */
                          }
/*                     printf(" retorno de la columna 22\n");       */
                     }
/*                 printf(" retorno de la columna 18\n");           */
                 }
                 else
                   readinline = 0;  
/*             printf(" retorno de la columna 14\n");               */
             }
             
             printf("already in getatomcoor 2\n");

           break;

        case 4:
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "PBC", 3)!= 0)
           {   
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "PBC", 3)!= 0)
           {
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }

           if (line == NULL)
             {
             getatomcoor = 0;
             break;
             }

/* need to put a counter, to stop while loop after reaching 100 000 lines in both readings above */
 
           sscanf(line,"%s %lf %lf %lf %lf %lf %lf", s1, &cell[1], &cell[2], &cell[3], &cell[4], &cell[5], &cell[6]);

           /* blanck line skipper */
           while (blanckline(line) == 0)
              {
                fgets(line, maxline, InProcesFile);
              }
           /*  end of blanck line skipper  */

           while (readinline == 1)
             {
               fgets(line, maxline, InProcesFile);
               if (line == NULL)
                {
                  getatomcoor = 0;
                  break;
                }
               while (blanckline(line) == 0)
                {
                  fgets(line, maxline, InProcesFile);
                }

/*               printf(" aviso0 %s\n", line);     */
               if (line != NULL)
                 {
                   if ( countwords(line, 0,strlen(line)-1) <= 3)
                     {
                       readinline = 0;
                     }
                     else
                     {
/*                       printf(" aviso1 %d\n", readinline);     
                       tmp_dummy = 1;
                       strcpy(s2 , grepnword(line, 0, 1));
                      printf(" palabra1 %s\n", s2);           
                       strcpy(s3 , "c");
                       printf(" palabra2 %s\n", s3);           
                       strcpy(s4 , grepnword(line, 0, 2));
                       printf(" palabra3 %s\n", s4);           
                       strcpy(s5 , grepnword(line, 0, 3));
*/
                       sscanf(line,"%s %s %s %s", s2, s4, s5, s6);
                       strcpy(s3 , "c");
/*                       printf(" read %s %s %s %s\n", s2, s3, s4, s5);          */
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            if ( countwords(line, 0,strlen(line)-1) <= 3)
                              readinline = 0;
                            else
                               tmp_dummy = 1;
/*                              strcpy(s6 , grepnword(line, 0, 4)); */
/*                            printf(" palabra5 %s\n", s6);      */
                            if (readinline != 0)
                               {
/*                                 printf(" aviso s6a %d\n", readinline);   */
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     {
                                       tmp_dummy = 1;
                                     }
                                   else
                                     {
                                       readinline = 0;
                                     }
                               }
                               else
                                 {
                                   readinline = 0;
/*                                   printf(" aviso s6b\n");      */
                                 }
                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {  
                                 readinline = 0;
                               }
                               
                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {
                                 readinline = 0;
                               }
                                                          
/*                            printf(" aviso readinline = %d\n", readinline);       */
                            if (readinline != 0)
                              {
                                atomnumber++;
                                sscanf(s4, "%lf", &xatom[atomnumber-1]);
/*                                printf("%f\t", xatom[atomnumber-1]);   */
                                sscanf(s5, "%lf", &yatom[atomnumber-1]);
/*                                printf(" %f\t", yatom[atomnumber-1]);  */
                                sscanf(s6, "%lf", &zatom[atomnumber-1]);
/*                                printf("%f\n", zatom[atomnumber-1]);   */
                                
                                tmp_dummy = 0;
                                for (tmp_dummy = 0; tmp_dummy <= strlen(s2) -1; tmp_dummy++)
                                   atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy];                                
                                /*
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                */
                                tmp_dummy = 1;
/*                                printf(" nombre atomo ?? %d %s\n", atomnumber-1, atomname[atomnumber-1]);  */
                               
                                atom_c_s[atomnumber-1] = s3[0];
/*                                printf(" %c\n", atom_c_s[atomnumber-1]);                               
                                                              
                                printf("ojo %d\t %s\t %c\t", atomnumber-1,  atomname[atomnumber-1], atom_c_s[atomnumber-1]);
                                printf("%f\t %f\t %f; %s\n", xatom[atomnumber-1], yatom[atomnumber-1], zatom[atomnumber-1], line); 
*/
                              }
                          }
                          else
                          {
                            printf(" llega a aqui ???");            
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
/*                          printf(" retorno de la columna 27\n");  */
                          }
/*                     printf(" retorno de la columna 22\n");       */
                     }
/*                 printf(" retorno de la columna 18\n");           */
                 }
                 else
                   readinline = 0;  
/*             printf(" retorno de la columna 14\n");               */
             }
             
             printf("already in getatomcoor 2\n");

           break;

        case 5:
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "_cell_length_a", 14)!= 0)
           {
//             printf("%s", line);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           // printf("test reading _cell_length_a %s\n", line);
           sscanf(line,"%s %s", s1, s2);
           // printf("test reading s2 %s\n", s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           strcpy(s1a, grepnmword(s2, 0, count2));
           strcpy(s1a, grepnmword(s2, 0, count2));
//           sscanf(grepnmword(s2, 0, count2), "%s", s1a);
           // printf("test reading s2 %s s1a %s and determining count2 %d\n", s2, s1a, count2);
           sscanf(s1a, "%lf", &cell[1]);
           // temporal para este problema!!!
           // sscanf(s2, "%lf", &cell[1]);
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           strcpy(s1, grepnmword(s2, 0, count2));
           sscanf(s1, "%lf", &cell[2]);
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           strcpy(s1, grepnmword(s2, 0, count2));
           sscanf(s1, "%lf", &cell[3]);
           fgets(line, maxline, InProcesFile); 
           sscanf(line,"%s %s", s1, s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           strcpy(s1, grepnmword(s2, 0, count2));
           sscanf(s1, "%lf", &cell[4]);
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           strcpy(s1, grepnmword(s2, 0, count2));
           sscanf(s1, "%lf", &cell[5]);
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           strcpy(s1, grepnmword(s2, 0, count2));
           sscanf(s1, "%lf", &cell[6]);
           printf("cell parameters %f %f %f %f %f %f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);
           
           while (strncmp(s1, "_atom_site_", 11) != 0)
           {
             chtmp1 = fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           printf("1 atom_site %s\n", s1);
           chtmp1 = fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           printf("2 atom_site %s\n", s1);
           tmp_dummy = 1;
           if (strncmp(s1, "_atom_site_fract_x", 18)!= 0)
           {
           	 tmp_dummy = 2;
		   } 
		   printf(" number of atom label columns %d\n", tmp_dummy);      
		   
           while (strncmp(s1, "_atom_site_fract_z", 18)!= 0)
           {
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           while (strncmp(s1, "_atom_site_", 11) == 0)
           {
             chtmp1 = fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           while ( (chtmp1 != NULL) && (strncmp(s1, "loop_", 5) != 0) )
           {
             if (strlen(line)>2)
             {
             	if (tmp_dummy == 2)
             	{
                  sscanf(line,"%s %s %lf %lf %lf %lf", s2, atomname[atomnumber], &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);
             	}
             	else
             	{
             	  sscanf(line,"%s %lf %lf %lf %lf", atomname[atomnumber], &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);  
             	}
                //printf("%d %s %.6f %.6f %.6f %.6f type input %d\n", atomnumber, atomname[atomnumber], xatom[atomnumber], yatom[atomnumber], zatom[atomnumber], atomcharge[atomnumber], tmp_dummy);
                atomnumber++;
             }
             chtmp1 = fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }

           strcpy(s3 , "c");
           for (count1 = 0; count1 < atomnumber; count1++)
           {
             atom_c_s[count1] = s3[0];
           }

           printf("already in getatomcoor 2\n");

           break;
        }
      
         printf("atomnumber %d\n", atomnumber); 
         printf("cell  %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]); 
         for (c1 = 0; c1 <= atomnumber-1; c1++)
               {
                 printf("%s %c %.5f %.5f %.5f\n", atomname[c1], atom_c_s[c1], xatom[c1], yatom[c1], zatom[c1]);  
               }

        fclose(InProcesFile);
        getatomcoor = 1;
        return getatomcoor;

}   /********************** end of getatomcoor ******************************/


/*****************************************************************************

this function gives the distances between two atoms in cartesian coordinates

*****************************************************************************/
double  atomsdistc(double xrotaxis1, double y1, double zrotaxis1, double x2, double y2, double z2)
{
    double distance;
        
    distance = sqrt((xrotaxis1-x2)*(xrotaxis1-x2) + (y1-y2)*(y1-y2) + (zrotaxis1-z2)*(zrotaxis1-z2));

    return distance;
}   /*************** end of atomscdist **********************/

/*****************************************************************************
             minor changes on the original function by Aileen Grey

this function gives the distances between two atoms of a periodic system
the function needs the lattice parameters given by an external variable
 extern double cell[N];

*****************************************************************************/
double  atomsdist(double x1, double y1, double z1, double x2, double y2, double z2)
{
        double atom[N];
        double o_atom[N];
        double per_atom[N];
        double ouratom[N];
        double o_ouratom[N];
        char line[201];
        char line2[201];
        char element[10];
        char type[10];
        double dist[10];
        double mag_dist;
        int i;
        int atom_num;

double  mag(double x, double y, double z);
double  coordtovec(int i, double cell[N], double per_atom[N]);
double  vectocoord(int i, double cell[N], double o_atom[N]);

 extern double cell[N];
        
        atom[1] = x1;
        atom[2] = y1;
        atom[3] = z1;
        ouratom[1] = x2;
        ouratom[2] = y2;
        ouratom[3] = z2;

        for (i=1;i<=3;i++)
          {
          if ((atom[i]-ouratom[i])>0.5)
             {
             per_atom[i]=(atom[i]-1.0);
             }
          else
             {
             if ((ouratom[i]-atom[i])>0.5)
               {
               per_atom[i]=(1.0+atom[i]);
               }
             else
               {
               per_atom[i]=atom[i];
               }
             }
          }
        for (i=1;i<=3;i++)
         {
         o_ouratom[i]=coordtovec(i, cell, ouratom);
         o_atom[i]=coordtovec(i, cell, per_atom);
         dist[i]=o_ouratom[i]-o_atom[i];
         }
         mag_dist = mag(dist[1], dist[2], dist[3]);

         return mag_dist;
}   /*************** end of the main body of atomsdist **********************/

/* function to find the magnitude of a vector */ 
double mag(double x, double y, double z)
{ 
          double mag;
          mag = sqrt((x*x)+(y*y)+(z*z)); 
          return mag;
}

/* function to calculate the cartesian coordinates relative to a */
/* orthonormal basis set with the x axis parallel to the x axis */
/* specified - units are in amstrongs */
double coordtovec(int i, double cell[N], double per_atom[N] )
{
  double a[N][N];
  double coordtovec;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));  
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3])); 

  coordtovec = (a[i][1]*per_atom[1])+(a[i][2]*per_atom[2])+(a[i][3]*per_atom[3]);

  return (coordtovec);
}

/*and function to change them back again */
double vectocoord(int i, double cell[N], double o_atom[N])
{
  double a[N][N];
  double acof[N][N];
  double ainv[N][N];
  double deta;
  double vectocoord;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3]));

  acof[1][1] = (a[2][2]*a[3][3])-(a[2][3]*a[3][2]);
  acof[1][2] = (a[1][3]*a[3][2])-(a[1][2]*a[3][3]);
  acof[1][3] = (a[1][2]*a[2][3])-(a[1][3]*a[2][2]);
  acof[2][1] = (a[2][3]*a[3][1])-(a[2][1]*a[3][3]);
  acof[2][2] = (a[1][1]*a[3][3])-(a[1][3]*a[3][1]);
  acof[2][3] = (a[1][3]*a[2][1])-(a[1][1]*a[2][3]);
  acof[3][1] = (a[2][1]*a[3][2])-(a[2][2]*a[3][1]);
  acof[3][2] = (a[1][2]*a[3][1])-(a[1][1]*a[3][2]);
  acof[3][3] = (a[1][1]*a[2][2])-(a[1][2]*a[2][1]);

  deta = (a[1][1]*acof[1][1])+(a[1][2]*acof[2][1])+(a[1][3]*acof[3][1]);

  ainv[1][1] = (1/deta)*(acof[1][1]);
  ainv[1][2] = (1/deta)*(acof[1][2]);
  ainv[1][3] = (1/deta)*(acof[1][3]);
  ainv[2][1] = (1/deta)*(acof[2][1]);
  ainv[2][2] = (1/deta)*(acof[2][2]);
  ainv[2][3] = (1/deta)*(acof[2][3]);
  ainv[3][1] = (1/deta)*(acof[3][1]);
  ainv[3][2] = (1/deta)*(acof[3][2]);
  ainv[3][3] = (1/deta)*(acof[3][3]);

  vectocoord = (ainv[i][1]*o_atom[1])+(ainv[i][2]*o_atom[2])+(ainv[i][3]*o_atom[3]);


  return (vectocoord);

}
/************************** end of atomsdist ********************************/

/*****************************************************************************
function gives 0 if the line not countain words, 1 else
*****************************************************************************/
int    blanckline(char line[Maxline])
{
        int    blanckline;
        int    c1, c2;
        
        blanckline = 0;
/*        printf(" in blancklin0 input:%s\n", line);     */
        if (strlen(line) == 1)
        {
/*           printf(" in blancklin1\n");                 */
           if (isspace(line[0]) == 0)
           {
             blanckline = 1;
           }
           else
           {
             blanckline = 0;
           }
        }
        else
        {
          for (c1 = 0; c1 <= strlen(line)-1;c1++)
            {
               if (isspace(line[c1]) == 0)
                 {
                   blanckline = 1;
                   c1 = strlen(line)+1;
                 }
            }       
        }


        return blanckline;

}    /***********************  end of blanckline  ****************************/

/*****************************************************************************
 function to calculate the cartesian coordinates relative to a 
 orthonormal basis set with the x axis parallel to the x axis 
 specified - units are in amstrongs 
 i = 1, 2 and 3 stand for  x, y and z, respectively
 
 original by Aileen Grey
*****************************************************************************/ 
double Cryst2Cartes(int i, double cell[N], double per_atom[N] )
{
  double a[N][N];
  double coordtovec;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));  
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3])); 

  coordtovec = (a[i][1]*per_atom[1])+(a[i][2]*per_atom[2])+(a[i][3]*per_atom[3]);

  return (coordtovec);
}  /**************  end of Cryst2Cartes ********************/ 

/*****************************************************************************
 function to calculate the crystallographic coordinates from cartesians
 i = 1, 2 and 3 stand for  x, y and z, respectively

 original by Aileen Grey
*****************************************************************************/
double Cartes2Cryst(int i, double cell[N], double per_atom[N] )
{
  double a[N][N];
  double acof[N][N];
  double ainv[N][N];
  double deta;
  double vectocoord;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3]));

  acof[1][1] = (a[2][2]*a[3][3])-(a[2][3]*a[3][2]);
  acof[1][2] = (a[1][3]*a[3][2])-(a[1][2]*a[3][3]);
  acof[1][3] = (a[1][2]*a[2][3])-(a[1][3]*a[2][2]);
  acof[2][1] = (a[2][3]*a[3][1])-(a[2][1]*a[3][3]);
  acof[2][2] = (a[1][1]*a[3][3])-(a[1][3]*a[3][1]);
  acof[2][3] = (a[1][3]*a[2][1])-(a[1][1]*a[2][3]);
  acof[3][1] = (a[2][1]*a[3][2])-(a[2][2]*a[3][1]);
  acof[3][2] = (a[1][2]*a[3][1])-(a[1][1]*a[3][2]);
  acof[3][3] = (a[1][1]*a[2][2])-(a[1][2]*a[2][1]);

  deta = (a[1][1]*acof[1][1])+(a[1][2]*acof[2][1])+(a[1][3]*acof[3][1]);

  ainv[1][1] = (1/deta)*(acof[1][1]);
  ainv[1][2] = (1/deta)*(acof[1][2]);
  ainv[1][3] = (1/deta)*(acof[1][3]);
  ainv[2][1] = (1/deta)*(acof[2][1]);
  ainv[2][2] = (1/deta)*(acof[2][2]);
  ainv[2][3] = (1/deta)*(acof[2][3]);
  ainv[3][1] = (1/deta)*(acof[3][1]);
  ainv[3][2] = (1/deta)*(acof[3][2]);
  ainv[3][3] = (1/deta)*(acof[3][3]);

  vectocoord = (ainv[i][1]*per_atom[1])+(ainv[i][2]*per_atom[2])+(ainv[i][3]*per_atom[3]);


  return (vectocoord);

}  /***************** end of Cartes2Cryst ************************/    


/****************************************************************************

simple function to gives absolute value of x

****************************************************************************/
double   absx(double x)
{
  double x1;
  
  if (x >= 0)
  {
    x1 = x;
  }
  else
  {
    x1 = -x;
  }
 
  return x1;
  
}  /***************** end of absx ************************/ 



/*****************************************************************************
this function gives the angles between three atoms of a periodic system
the function needs the lattice parameters given by an external variable
extern double cell[N];

the formula is that by J. Buerger, appeared in his book "", chapter 23, pp. 629

atom1 is the centre of the angle

******************************************************************************/
double  atomsangle(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33)
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    double s12, s13;
    double deltax12, deltax13, deltay12, deltay13, deltaz12, deltaz13;
    double angle1;

    /**  first put all the atoms in the same unit cell, the next to the (0,0,0) one at the positive site **/
    if ( (x11 > 1.0) || (x11 < 0.0) )
    {
      if (x11 > 1.0)
      {
       x1 = x11 - 1.0;
      }
      else
      {
        x1 = x11 + 1.0;
      } 
    }
    else
    {
     x1 = x11;
    }
    
    if ( (x22 > 1.0) || (x22 < 0.0) )
    {
      if (x22 > 1.0)
      {
       x2 = x22 - 1.0;
      }
      else
      {
        x2 = x22 + 1.0;
      } 
    }
    else
    {
     x2 = x22;
    }
    
    if ( (x33 > 1.0) || (x33 < 0.0) )
    {
      if (x33 > 1.0)
      {
       x3 = x33 - 1.0;
      }
      else
      {
        x3 = x33 + 1.0;
      } 
    }
    else
    {
     x3 = x33;
    }
    
    if ( (y11 > 1.0) || (y11 < 0.0) )
    {
      if (y11 > 1.0)
      {
       y1 = y11 - 1.0;
      }
      else
      {
        y1 = y11 + 1.0;
      } 
    }
    else
    {
     y1 = y11;
    }
    
    if ( (y22 > 1.0) || (y22 < 0.0) )
    {
      if (y22 > 1.0)
      {
       y2 = y22 - 1.0;
      }
      else
      {
        y2 = y22 + 1.0;
      } 
    }
    else
    {
     y2 = y22;
    }
    
    if ( (y33 > 1.0) || (y33 < 0.0) )
    {
      if (y33 > 1.0)
      {
       y3 = y33 - 1.0;
      }
      else
      {
        y3 = y33 + 1.0;
      } 
    }
    else
    {
     y3 = y33;
    }
    
    if ( (z11 > 1.0) || (z11 < 0.0) )
    {
      if (z11 > 1.0)
      {
       z1 = z11 - 1.0;
      }
      else
      {
        z1 = z11 + 1.0;
      } 
    }
    else
    {
     z1 = z11;
    }
    
    if ( (z22 > 1.0) || (z22 < 0.0) )
    {
      if (z22 > 1.0)
      {
       z2 = z22 - 1.0;
      }
      else
      {
        z2 = z22 + 1.0;
      } 
    }
    else
    {
     z2 = z22;
    }
    
    if ( (z33 > 1.0) || (z33 < 0.0) )
    {
      if (z33 > 1.0)
      {
       z3 = z33 - 1.0;
      }
      else
      {
        z3 = z33 + 1.0;
      } 
    }
    else
    {
     z3 = z33;
    }
    /**  end of  put all the atoms in the same unit cell, the next to the (0,0,0) one at the positive site **/
    
    /**  first put the second atom coordinates as the closer distances of the first atom **/
    if ( (x1 - x2 > 0.5) || (x1 - x2 < -0.5) )
    {
      if (x1 - x2 > 0.5) 
      {
        x2 = x2 + 1.0;
      }
      else
      {
        x2 = x2 - 1.0;
      }
    }
    
    if ( (y1 - y2 > 0.5) || (y1 - y2 < -0.5) )
    {
      if (y1 - y2 > 0.5)
      {
        y2 = y2 + 1.0;
      }
      else
      {
        y2 = y2 - 1.0;
      }
    }
    
    if ( (z1 - z2 > 0.5) || (z1 - z2 < -0.5) )
    {
      if (z1 - z2 > 0.5)
      {
        z2 = z2 + 1.0;
      }
      else
      {
        z2 = z2 - 1.0;
      }
    }
    /**  end of put the second atom coordinates as the closer distances of the first atom **/
    
    /**  first put the third atom coordinates as the closer distances of the first atom **/
    if ( (x1 - x3 > 0.5) || (x1 - x3 < -0.5) )
    {
      if (x1 - x3 > 0.5) 
      {
        x3 = x3 + 1.0;
      }
      else
      {
        x3 = x3 - 1.0;
      }
    }
    
    if ( (y1 - y3 > 0.5) || (y1 - y3 < -0.5) )
    {
      if (y1 - y3 > 0.5)
      {
        y3 = y3 + 1.0;
      }
      else
      {
        y3 = y3 - 1.0;
      }
    }
    
    if ( (z1 - z3 > 0.5) || (z1 - z3 < -0.5) )
    {
      if (z1 - z3 > 0.5)
      {
        z3 = z3 + 1.0;
      }
      else
      {
        z3 = z3 - 1.0;
      }
    }

    /**  end of put the third atom coordinates as the closer distances of the first atom **/
    
    /**     this is the key segment of the function       **/
    s12 = atomsdist(x1, y1, z1, x2, y2, z2);
    s13 = atomsdist(x1, y1, z1, x3, y3, z3);    
    
/*    printf(" distance 1 %f    distance 2 %f\n", s12, s13); */
    
    deltax12 = x1 - x2;
    deltax13 = x1 - x3;
    deltay12 = y1 - y2;
    deltay13 = y1 - y3;
    deltaz12 = z1 - z2;
    deltaz13 = z1 - z3;
/*    
    printf(" deltax12 %f    deltax13  %f\n", deltax12, deltax13);
    printf(" deltay12 %f    deltay13  %f\n", deltay12, deltay13);
    printf(" deltaz12 %f    deltaz13  %f\n", deltaz12, deltaz13);
*/    
    angle1 = deltax12*deltax13*cell[1]*cell[1] + deltay12*deltay13*cell[2]*cell[2] + deltaz12*deltaz13*cell[3]*cell[3];
/*    printf(" 1            angle %f\n", angle1); */
    angle1 = angle1 + (deltax12*deltay13 + deltay12*deltax13)*cell[1]*cell[2]*cos(cell[6]*M_PI/180.0);
/*    printf(" 2            angle %f\n", angle1); */
    angle1 = angle1 + (deltaz12*deltax13 + deltax12*deltaz13)*cell[1]*cell[3]*cos(cell[5]*M_PI/180.0);
/*    printf(" 3            angle %f\n", angle1); */
    angle1 = angle1 + (deltay12*deltaz13 + deltaz12*deltay13)*cell[2]*cell[3]*cos(cell[4]*M_PI/180.0);
/*    printf(" 4            angle %f\n", angle1); */
    
    angle1 = angle1/(s12*s13);
    if (angle1 > 1.0) 
     { 
       angle1 = 1.0 ; 
     }
    else 
     { if ( angle1 < -1.0 ) 
       { angle1 = -1.0;
       }
     }
/*    printf(" cosin of the angle %f\n", angle1); */
    
    angle1 = acos(angle1)*180.0/M_PI;
/*    printf("              angle %f\n", angle1); */
    
    
    return angle1;
}
/************************** end of atomsangle ********************************/


/*****************************************************************************
this function gives the angles between three atoms of a periodic system
the function needs the lattice parameters given by an external variable
extern double cell[N];

the formula is that by J. Buerger, appeared in his book "", chapter 23, pp. 629

atom1 is the centre of the angle

******************************************************************************/
double  atomsangle_test(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33)
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    double s12, s13;
    double deltax12, deltax13, deltay12, deltay13, deltaz12, deltaz13;
    double angle1;


    printf(" atomsangle_test ms1 input coord atom1 %f %f %f atom2 %f %f %f atom3 %f %f %f\n", x11, y11, z11, x22, y22, z22, x33, y33, z33);

    /**  first put all the atoms in the same unit cell, the next to the (0,0,0) one at the positive site **/
    if ( (x11 > 1.0) || (x11 < 0.0) )
    {
      if (x11 > 1.0)
      {
       x1 = x11 - 1.0;
      }
      else
      {
        x1 = x11 + 1.0;
      } 
    }
    else
    {
     x1 = x11;
    }
    
    if ( (x22 > 1.0) || (x22 < 0.0) )
    {
      if (x22 > 1.0)
      {
       x2 = x22 - 1.0;
      }
      else
      {
        x2 = x22 + 1.0;
      } 
    }
    else
    {
     x2 = x22;
    }
    
    if ( (x33 > 1.0) || (x33 < 0.0) )
    {
      if (x33 > 1.0)
      {
       x3 = x33 - 1.0;
      }
      else
      {
        x3 = x33 + 1.0;
      } 
    }
    else
    {
     x3 = x33;
    }
    
    if ( (y11 > 1.0) || (y11 < 0.0) )
    {
      if (y11 > 1.0)
      {
       y1 = y11 - 1.0;
      }
      else
      {
        y1 = y11 + 1.0;
      } 
    }
    else
    {
     y1 = y11;
    }
    
    if ( (y22 > 1.0) || (y22 < 0.0) )
    {
      if (y22 > 1.0)
      {
       y2 = y22 - 1.0;
      }
      else
      {
        y2 = y22 + 1.0;
      } 
    }
    else
    {
     y2 = y22;
    }
    
    if ( (y33 > 1.0) || (y33 < 0.0) )
    {
      if (y33 > 1.0)
      {
       y3 = y33 - 1.0;
      }
      else
      {
        y3 = y33 + 1.0;
      } 
    }
    else
    {
     y3 = y33;
    }
    
    if ( (z11 > 1.0) || (z11 < 0.0) )
    {
      if (z11 > 1.0)
      {
       z1 = z11 - 1.0;
      }
      else
      {
        z1 = z11 + 1.0;
      } 
    }
    else
    {
     z1 = z11;
    }
    
    if ( (z22 > 1.0) || (z22 < 0.0) )
    {
      if (z22 > 1.0)
      {
       z2 = z22 - 1.0;
      }
      else
      {
        z2 = z22 + 1.0;
      } 
    }
    else
    {
     z2 = z22;
    }
    
    if ( (z33 > 1.0) || (z33 < 0.0) )
    {
      if (z33 > 1.0)
      {
       z3 = z33 - 1.0;
      }
      else
      {
        z3 = z33 + 1.0;
      } 
    }
    else
    {
     z3 = z33;
    }
    /**  end of  put all the atoms in the same unit cell, the next to the (0,0,0) one at the positive site **/
    
    /**  first put the second atom coordinates as the closer distances of the first atom **/
    if ( (x1 - x2 > 0.5) || (x1 - x2 < -0.5) )
    {
      if (x1 - x2 > 0.5) 
      {
        x2 = x2 + 1.0;
      }
      else
      {
        x2 = x2 - 1.0;
      }
    }
    
    if ( (y1 - y2 > 0.5) || (y1 - y2 < -0.5) )
    {
      if (y1 - y2 > 0.5)
      {
        y2 = y2 + 1.0;
      }
      else
      {
        y2 = y2 - 1.0;
      }
    }
    
    if ( (z1 - z2 > 0.5) || (z1 - z2 < -0.5) )
    {
      if (z1 - z2 > 0.5)
      {
        z2 = z2 + 1.0;
      }
      else
      {
        z2 = z2 - 1.0;
      }
    }
    /**  end of put the second atom coordinates as the closer distances of the first atom **/
    
    /**  first put the third atom coordinates as the closer distances of the first atom **/
    if ( (x1 - x3 > 0.5) || (x1 - x3 < -0.5) )
    {
      if (x1 - x3 > 0.5) 
      {
        x3 = x3 + 1.0;
      }
      else
      {
        x3 = x3 - 1.0;
      }
    }
    
    if ( (y1 - y3 > 0.5) || (y1 - y3 < -0.5) )
    {
      if (y1 - y3 > 0.5)
      {
        y3 = y3 + 1.0;
      }
      else
      {
        y3 = y3 - 1.0;
      }
    }
    
    if ( (z1 - z3 > 0.5) || (z1 - z3 < -0.5) )
    {
      if (z1 - z3 > 0.5)
      {
        z3 = z3 + 1.0;
      }
      else
      {
        z3 = z3 - 1.0;
      }
    }

    /**  end of put the third atom coordinates as the closer distances of the first atom **/

    printf(" atomsangle_test ms2 trans coord atom1 %f %f %f atom2 %f %f %f atom3 %f %f %f\n", x11, y11, z11, x22, y22, z22, x33, y33, z33);
    
    /**     this is the key segment of the function       **/
    s12 = atomsdist(x1, y1, z1, x2, y2, z2);
    s13 = atomsdist(x1, y1, z1, x3, y3, z3);    
    
    printf(" atomsangle_test ms3 dist12 %f dist13 %f\n", s12, s13);

/*    printf(" distance 1 %f    distance 2 %f\n", s12, s13); */
    
    deltax12 = x1 - x2;
    deltax13 = x1 - x3;
    deltay12 = y1 - y2;
    deltay13 = y1 - y3;
    deltaz12 = z1 - z2;
    deltaz13 = z1 - z3;
/*    
    printf(" deltax12 %f    deltax13  %f\n", deltax12, deltax13);
    printf(" deltay12 %f    deltay13  %f\n", deltay12, deltay13);
    printf(" deltaz12 %f    deltaz13  %f\n", deltaz12, deltaz13);
*/    
    angle1 = deltax12*deltax13*cell[1]*cell[1] + deltay12*deltay13*cell[2]*cell[2] + deltaz12*deltaz13*cell[3]*cell[3];
    printf(" atomsangle_test ms4 tmp1 angle1 %f\n", angle1);
/*    printf(" 1            angle %f\n", angle1); */
    angle1 = angle1 + (deltax12*deltay13 + deltay12*deltax13)*cell[1]*cell[2]*cos(cell[6]*M_PI/180.0);
    printf(" atomsangle_test ms4 tmp2 angle1 %f\n", angle1);
/*    printf(" 2            angle %f\n", angle1); */
    angle1 = angle1 + (deltaz12*deltax13 + deltax12*deltaz13)*cell[1]*cell[3]*cos(cell[5]*M_PI/180.0);
    printf(" atomsangle_test ms4 tmp3 angle1 %f\n", angle1);
/*    printf(" 3            angle %f\n", angle1); */
    angle1 = angle1 + (deltay12*deltaz13 + deltaz12*deltay13)*cell[2]*cell[3]*cos(cell[4]*M_PI/180.0);
    printf(" atomsangle_test ms4 tmp4 angle1 %f\n", angle1);
/*    printf(" 4            angle %f\n", angle1); */
    
    angle1 = angle1/(s12*s13);
    printf(" atomsangle_test ms4 tmp5 angle1 %f\n", angle1);
    if (angle1 > 1.0) 
    { 
      angle1 = 1.0 ; 
      printf(" atomsangle_test ms4 tmp6+ OJO angle1 %f\n", angle1);
    }
    else 
    { 
      if ( angle1 < -1.0 ) 
      { 
        angle1 = -1.0;
        printf(" atomsangle_test ms4 tmp6- OJO angle1 %f\n", angle1);  
      }
    }
/*    printf(" cosin of the angle %f\n", angle1); */
    
    angle1 = acos(angle1)*180.0/M_PI;
    printf(" atomsangle_test ms4 FINAL angle1 %f\n", angle1);
/*    printf("              angle %f\n", angle1); */
    
    
    return angle1;
}
/************************** end of atomsangle ********************************/

/*****************************************************************************
this function gives the angles between three atoms in a cartesian non-periodic
system

******************************************************************************/
double  atomsanglec(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33)
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    double s12, s13, r12xr13;
    double deltax12, deltax13, deltay12, deltay13, deltaz12, deltaz13;
    double angle1;

    s12 = atomsdistc(x11, y11, z11, x22, y22, z22);
    s13 = atomsdistc(x11, y11, z11, x33, y33, z33);    
    
/*    printf(" distance 1 %f    distance 2 %f\n", s12, s13); */
    
    r12xr13 = (x22-x11)*(x33-x11) + (y22-y11)*(y33-y11) + (z22-z11)*(z33-z11);

    angle1 = r12xr13/(s12*s13);
/*    printf(" cosin of the angle %f\n", angle1); */
    
    angle1 = acos(angle1)*180.0/M_PI;
/*    printf("              angle %f\n", angle1); */
    
    
    return angle1;
}
/************************** end of atomsanglec ********************************/

/****************************************************************************
  this function gives the nth word in a string
*****************************************************************************/
char  *grepnword(char *line1, int position_line, int nthword)
{
        int    c1, c2, c3, tempcount;
        char  grepnword1[Maxline];

        c2 = position_line;
        
        for (tempcount = 1; tempcount <= nthword; tempcount++)
        {
          c1 = c2;
          while ( (line1[c1] == ' ' )  && (c1 <= (strlen(line1)-1) ) )
             {
             c1++;
             }
          c2 = c1+1;
          while ( (line1[c2] != ' ')  && (c2 <= (strlen(line1)-1) ) )
             {
               c2++;
             }
          if (c2 == (strlen(line1)))
             {
               if (tempcount < nthword)
               {
                 tempcount = nthword+2;
                 strcpy(grepnword1, "");
               }
             }
        }

        if (tempcount != nthword+2)
        {
          for (c3 = 0; c3 < c1; c3++)
            line1++;
          sscanf(line1,"%s", grepnword1);
        }

        return grepnword1;

} /***********************  end of grepnword *******************************/


/*****************************************************************************
function to grep a string (line2) in a file,
returns 0  if the string is in the file, and
returns -1 if the string isn't in the file  
also returns as 'line' the line read, if exits !
*****************************************************************************/
int    strgrepf(char line1[400], char line2[400], int maxline1, FILE *filename1) 
{   
    char  *chtmp1;
    int   auxcount1;
    int   strgrepf;
    
    auxcount1 = 0;
    chtmp1 = fgets(line1, maxline1, filename1);
    while ((strstr(line1, line2) == NULL) && (chtmp1 != NULL))
        {                 
          chtmp1 = fgets(line1, maxline1, filename1);
        }
    if (chtmp1 != NULL)    
       strgrepf = 0;
    else 
       strgrepf = -1;

   return strgrepf;            
}    /***********************  end of strgrepf  ****************************/


/****************************************************************************
     this function gives the number of words in a line between
     two given points; count1pos, count2pos.
*****************************************************************************/
int    countwords(char line1[Maxline], int count1pos, int count2pos)
{
        int   c1, state_word;
        int   countwords;

        
        if ( isspace(line1[count1pos]) == 0) 
          {
            state_word = 1;
            countwords = 1;
          }
          else
          {
            state_word = 0;
            countwords = 0;
          }
          
          c1 = count1pos-1;
          while (c1 < count2pos)
            {
              c1++;
              if (isspace(line1[c1]) == 0)
                {
                  if (state_word == 0)
                    {
                      state_word = 1;
                      countwords++;
                    }
                }
                else
                {
                  if (state_word != 0)
                    state_word = 0;
                }   
            }

        return countwords;

}   /********************** end of countwords ******************************/

/****************************************************************************
  this function gives the word located between position n and m of a string
*****************************************************************************/
char  *grepnmword(char *line1, int position_ini, int position_final)
{
        int    c1, c2, c3, tempcont;
        char  grepnmword1[Maxline];
        char  grepnmword2[Maxline];
        char  grepnmword3[Maxline];
        char  ch1;
        int   ch_int1, ch_int2;

        // printf("test grepnmword reading line1 %s position_ini %d position_final %d\n", line1, position_ini, position_final);

        strcpy(grepnmword2, "");


        if ((position_ini > strlen(line1)) || (position_ini > position_final) )
        {
          strcpy(grepnmword2, "");
        }
        else
        {
          if (position_final > strlen(line1))
            position_final = strlen(line1);
          for (c3 = 0; c3 < position_ini; c3++)
            line1++;

          ch_int1 = line1[0];
          sprintf(grepnmword2, "%c", ch_int1);
          sprintf(grepnmword3, "%c", ch_int1);

          c2 = position_final - position_ini;
          for (c1 = 1; c1 < c2; c1++)
          {
            line1++;
            ch_int1 = line1[0];
            sprintf(grepnmword3, "%c", ch_int1);
            strcat(grepnmword2, grepnmword3);
          }
        }

        // printf("test grepnmword reading grepnmword2 %s\n", grepnmword2);

        return grepnmword2;

} /***********************  end of grepnmword *******************************/
                                                                               
/****************************************************************************
     this function uses getatomcoor to extract the coordinates and to create
     a list of the singleT (stored in singleTatoms)
*****************************************************************************/
void   FindSingleT(char CoordFNamex[Maxline], int CoorFTypex)
{
    int    count1, count2, count3, count4;
    int    tmpAtoms[N];
    double disttmp;
    char   s1[Maxline], s2[Maxline];

    printf(" ms8 %s %d\n", CoordFNamex, CoorFTypex);
    
    countsingleT = -1; countoxygen = -1;
//    count3 = (MaxAtomNum / 3) + 1;
    for (count1 = 0; count1 < MaxAtomNum2; count1++)
    {
      for (count2 = 0; count2 < 5; count2++)
      {
        singleTatoms[count2][count1] = -1;
      }
    }
    for (count1 = 0; count1 < MaxAtomNum3; count1++)
    {
      for (count2 = 0; count2 < 3; count2++)
      {
        oxyg_neighb[count2][count1] = -1;
      }
    }
    
    printf(" already in FindSingleT 1, ready to process file %s of type %d\n", CoordFNamex, CoorFTypex);
    
    count1 = getatomcoor(CoordFNamex, CoorFTypex);
    
    //count2 = atomnumber;
    //count3 = atomnumber;
    //for (count1 = 0; count1 < count2; count1++)
    //{
    //  strcpy(s1, atomname[count1]);
    //  if ( (s1[0] !='s') && (s1[0] !='S') && (s1[0] !='A') && (s1[0] !='A') && (s1[0] !='G') && (s1[0] !='g') && (s1[0] !='o') && (s1[0] !='O') ) 
    //  {
    //    atomnumber = count1;
    //    count1     = count2 + 10;
    //  }
    //}
    //
    //if (count3 != atomnumber)
    //{
    //  printf("         in FindSingleT, atomnumber change old %d new %d\n", count3, atomnumber);
    //}
    //else
    //{
    //  printf("         in FindSingleT, atomnumber not change\n");
    //}
    
    printf(" cell parameters %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]); 
    
/* 
    printf(" atom 0 is %s %4.5f %4.5f %4.5f \n", atomname[0], xatom[0], yatom[0], zatom[0]);
    printf(" atom 1 is %s %4.5f %4.5f %4.5f \n", atomname[1], xatom[1], yatom[1], zatom[1]);
    printf(" atom 2 is %s %4.5f %4.5f %4.5f \n", atomname[2], xatom[2], yatom[2], zatom[2]);
    printf(" atom 3 is %s %4.5f %4.5f %4.5f \n", atomname[3], xatom[3], yatom[3], zatom[3]);
    printf(" atom 4 is %s %4.5f %4.5f %4.5f \n", atomname[4], xatom[4], yatom[4], zatom[4]);
    printf(" atomnumber %d \n", atomnumber);
*/    
    printf(" already in FindSingleT 2, getatomcoor done!\n");

    for (count1 = 0; count1 < atomnumber; count1++)
    {
      tmpAtoms[0] = count1; 
      tmpAtoms[1] = -1;    tmpAtoms[2] = -1;
      tmpAtoms[3] = -1;    tmpAtoms[4] = -1;

      strcpy(s1, atomname[count1]);
      if ( (s1[0] =='s') || (s1[0] =='S') || (s1[0] =='A') || (s1[0] =='A') || (s1[0] =='G') || (s1[0] =='g') || (s1[0] =='o') || (s1[0] =='O') )    
	  {   
        count3 = 0; 
        for (count2 = 0; count2 < atomnumber; count2++)
        {
          strcpy(s2, atomname[count2]);
          if ( (s2[0] =='s') || (s2[0] =='S') || (s2[0] =='A') || (s2[0] =='A') || (s2[0] =='G') || (s2[0] =='g') || (s2[0] =='o') || (s2[0] =='O') )    
	      {
            disttmp = atomsdist(xatom[count1], yatom[count1], zatom[count1], xatom[count2], yatom[count2], zatom[count2]);
            if ((disttmp<=TOdistmax) && (disttmp>=TOdistmin))
            {
              count3++;
              tmpAtoms[count3] = count2;
            }
          }
        }
        printf(" current atom %d with %d bonds\n", count1, count3);
        if (count3 == 4)
        {
          countsingleT++;
          singleTatoms[0][countsingleT] = tmpAtoms[0];
          singleTatoms[1][countsingleT] = tmpAtoms[1];
          singleTatoms[2][countsingleT] = tmpAtoms[2];
          singleTatoms[3][countsingleT] = tmpAtoms[3];
          singleTatoms[4][countsingleT] = tmpAtoms[4];
        
          printf(" atom 0 of singleTatoms %d is %d\n", countsingleT, singleTatoms[0][countsingleT]);
          printf(" atom 1 of singleTatoms %d is %d\n", countsingleT, singleTatoms[1][countsingleT]);
          printf(" atom 2 of singleTatoms %d is %d\n", countsingleT, singleTatoms[2][countsingleT]);
          printf(" atom 3 of singleTatoms %d is %d\n", countsingleT, singleTatoms[3][countsingleT]);
          printf(" atom 4 of singleTatoms %d is %d\n", countsingleT, singleTatoms[4][countsingleT]);
        }

        if (count3 == 2)
        {
          countoxygen++;
          oxyg_neighb[0][countoxygen] = tmpAtoms[0];
          oxyg_neighb[1][countoxygen] = tmpAtoms[1];
          oxyg_neighb[2][countoxygen] = tmpAtoms[2];
        
          printf(" atom 0 of oxyg_neighb %d is %d\n", countoxygen, oxyg_neighb[0][countoxygen]);
          printf(" atom 1 of oxyg_neighb %d is %d\n", countoxygen, oxyg_neighb[1][countoxygen]);
          printf(" atom 2 of oxyg_neighb %d is %d\n", countoxygen, oxyg_neighb[2][countoxygen]);
        }
      }
    }
 
    printf("\n");
 
}   /********************** end of FindSingleT ******************************/

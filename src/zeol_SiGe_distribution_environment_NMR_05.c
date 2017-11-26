/****************************************************************************

-----------------------------------------------------------------------------
----------------    zeol_SiGe_distribution_environment_NMR_01   -------------
-----------        A.Rabdel Ruiz-Salvador(a)                    -------------
----           (a) Univ. Pablo de Olabvide, Seville, Spain               ----
--------------------------    June 2017         -----------------------------                              
-----------------------------------------------------------------------------

  The codes takes the pure silica/germania zeolite framework and put 
  heteroatoms in the place existing in reference structures
  Files names are listed in a file that is read by the code
    
  crash_label is an error message flag: = 0 if error ocurrs, > 0 else 
*****************************************************************************/

# include <math.h>
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>

//# define MaxAtomNum   500000
# define MaxAtomNum   1000
# define MaxAtomNum1  200000
# define MaxAtomNum2  350000
# define MaxAtomNum3  1000
# define MaxAtomNum4  1000
# define MaxAtomNum5  50

# define Maxline      300
# define N            20  
# define TOidealDist  1.605
# define TidealAngle  109.47122
  
FILE  *InpF;
FILE  *OutpF;
char   main_Output_FName[Maxline];
double xx1[MaxAtomNum], yy1[MaxAtomNum], zz1[MaxAtomNum], xatom[MaxAtomNum], yatom[MaxAtomNum], zatom[MaxAtomNum];
char   atomname[MaxAtomNum][8];
char   atom_c_s[MaxAtomNum];
int    CoorFType_monoatom_input, CoorFType_heteroatom_input[MaxAtomNum4];
long int atomnumber, atomnumber1;
double cell[N], cellRef[N];
char   monoatom_input_FName[Maxline], heteroatom_input_FName[MaxAtomNum4][Maxline], main_files_nameF[Maxline];
//double tolerance_dist;
double atomcharge[MaxAtomNum];
//long int site_atom[MaxAtomNum], printing_atom[MaxAtomNum], atom_input_order[MaxAtomNum];
int    number_input_files;
char   gulp_topFN[Maxline], gulp_bottomFN[Maxline];
int    number_Ge, Ge_atom[Maxline][N]; // 1 D4R (1-6); 2 atoms 1-8 in D4R 
double energy_config[MaxAtomNum4];
int    D4Ratoms[N][N][N]; // 1 D4R (1-6); 2 atoms 0=F 1-8 T(Si=1,Ge=2); 3 0-8 0=F, 1=itself, 2-4 neighbour, 5-7 diagonal face, 8 main diagonal 
int    Ge_in_D4R[N][N];   // 1 D4R (1-6); 2 0=number of Ge


void   Initialization();
void   MainJob();
int    getatomcoor(char coordFName[Maxline], int CoorFType);
char  *grepnword(char *line1, int position_line, int nthword);
char  *grepnmword(char *line1, int position_ini, int position_final);
int    strgrepf(char line1[], char line2[],int maxline1, FILE *filename1);
int    countwords(char line1[Maxline], int count1pos, int count2pos);
double absx(double x);
double atomsdist(double x1, double y1, double z1, double x2, double y2, double z2);

int main(argc, argv)
        int argc;
        char **argv;
{
int   count1, count2, count3, crash_label; 

     crash_label = 0.0;
     InpF = fopen(argv[1],"r");
     
     Initialization();
     
     printf(" passed over Initialization\n");

     printf("\n");
     
     //printf(" SingleTFName  %s  CoorFType  %d\n", SingleTFName, CoorFType);
     
     MainJob();
 
     printf(" thanks for using zeol_SiGe_distribution_environment_NMR_05 \n");
    
}     /******************************  end of main  **************************/

/******************************************************************************
function to initiate the job; reading general data and initialising
general variables.
******************************************************************************/
void    Initialization()
{   
    int  maxline, count1;
    char line[Maxline], s1[Maxline], s2[Maxline];

    maxline = Maxline;     
    
    fgets(line, maxline, InpF);
    sscanf(line,"%s %d %d", main_files_nameF, &number_input_files, &number_Ge);
    // name of the file that contains the files name: input heteroatom, file type  and output heteroatom
    
    fgets(line, maxline, InpF);
    sscanf(line,"%s %d", monoatom_input_FName, &CoorFType_monoatom_input);
    // name and type (gin, gout, xtl, car, cif) of the structure with bad order
	  
    fgets(line, maxline, InpF);
    sscanf(line,"%s", main_Output_FName);
    // Output (cif) good ordered atoms file nane 
        
    fclose(InpF);
    
    printf(" ms1 %s %d %d\n", main_files_nameF, number_input_files, number_Ge);
    printf(" ms2 %s %d\n", monoatom_input_FName,  CoorFType_monoatom_input);
    printf(" ms3 %s\n", main_Output_FName);
    
}    /***********    end of Initialization     ************/

/******************************************************************************
function to do the Main Job
******************************************************************************/
void   MainJob()
{
    char   line[Maxline];
    int    count3, count4, count5, count6, count7, maxline;
    double disttmp;
    int    count1, count2, atomnumber_good_input;
    char   tmp_inp_FN[Maxline], tmp_outp_FN[Maxline]; 
    char   s1[Maxline], s2[Maxline], s3[Maxline], currentFname[Maxline];
    char * pch;
    int    tmp_char, atomnumberRefi, Tatomnumber;
    int    tmp_D4R[N], count8;
    double tmp_dist[N], distmin, used_tmp_D4R[N];
    int    D4R_stat[30], MainCount, MainCount1, MainCount2, MainCount3;

    maxline = Maxline;

    count4 = getatomcoor(monoatom_input_FName, CoorFType_monoatom_input);
 
    for (count3 = 1; count3 <= 6; count3++)
    {
      cellRef[count3] = cell[count3];
    }

    // storing coordinates of the monoatom input structure
    for (count1 = 0; count1 < atomnumber; count1++)
    {
      xx1[count1] = xatom[count1];
      yy1[count1] = yatom[count1];
      zz1[count1] = zatom[count1];
    }
    // for checking
    printf("monoatom input structure atomnumber %d\n", atomnumber);
    for (count1 = 0; count1 < atomnumber; count1++)
    {
     printf("monoatom input structure: %s %c %f %f %f\n", atomname[count1], atom_c_s[count1], xx1[count1], yy1[count1], zz1[count1]);
    }

    // initializing D4Ratoms
    for (count1 = 0; count1 < N; count1++)
    {
      for (count2 = 0; count2 < N; count2++)
      { 
        for (count3 = 0; count3 < N; count3++)
        {
          D4Ratoms[count1][count2][count3] = -1;
        }
      }
    }

    // proccessing D4Ratoms
    for (count1 = atomnumber-6; count1 < atomnumber; count1++)
    {
      count5 = count1 - atomnumber + 7;
      D4Ratoms[count5][0][0] = count1;
      count6  = 0;
      count4  = -1;
      distmin = 10000.0;
      while ( (count6 != 8) && (count6 >= 0) )
      {
        for (count2 = 0; count2 < atomnumber; count2++)
        {
          disttmp = atomsdist(xx1[count1], yy1[count1], zz1[count1], xx1[count2], yy1[count2], zz1[count2]);
          //printf(" ms3a tmp_D4R %d at F %d vertex %d scanning T_atom %d disttmp %f\n", count5, count1, count6, count2, disttmp);
          if ( (disttmp <= 4.0) && (disttmp >= 0.01) )
          {
            count6++;
            tmp_D4R[count6]  = count2;
            if (disttmp < distmin)
            {
              distmin = disttmp;
              count4  = count2;
            }
            printf(" ms4a tmp_D4R %d at F %d vertex %d atom %d tmp_dist %f\n", count5, count1, count6, count2, disttmp);
            if (count6 == 8)
            {
              count2 = atomnumber + 1;
            }
          }
        }
      }
//    }
//    if (count6 > 10000) // OJO aqui esta el corte de test
//    {
      // finding D4Ratoms[count5][count7][1] 
      count7 = 1;
      D4Ratoms[count5][count7][1] = count4;
      for (count3 = 1; count3 <= 8; count3++)
      {
        if (count4 == tmp_D4R[count3])
        {
          used_tmp_D4R[count3] = 1;
        }
        else
        {
          used_tmp_D4R[count3] = 0;
        }
      }
      count7++;
      while (count7 < 9)
      {
        distmin = 10000.0;
        count4  = -1;
        for (count3 = 1; count3 <= 8; count3++)
        {
          if (used_tmp_D4R[count3] == 0)
          {
            count6 = tmp_D4R[count3];
            disttmp = atomsdist(xx1[count1], yy1[count1], zz1[count1], xx1[count6], yy1[count6], zz1[count6]);
            if (disttmp < distmin)
            {
              distmin = disttmp;
              count4  = count3;
            }
          }
        }
        used_tmp_D4R[count4] = 1;
	D4Ratoms[count5][count7][1] = tmp_D4R[count4];         
        count7++;
      }

      // finding D4Ratoms[count5][count7][2-8]
      for (count7 = 1; count7 <= 8; count7++)
      {
        count2 = D4Ratoms[count5][count7][1];
        for (count3 = 1; count3 <= 8; count3++)
        {
          if (count2 == tmp_D4R[count3])
          {
            used_tmp_D4R[count3] = 1;
          }
          else
          {
            used_tmp_D4R[count3] = 0;
          }
        } 
        for (count8 = 2; count8 <= 8; count8++)
        {
          distmin = 10000.0;
          count4  = -1;
          for (count3 = 1; count3 <= 8; count3++)
          {
            if (used_tmp_D4R[count3] == 0)
            {
              count6 = tmp_D4R[count3];
              disttmp = atomsdist(xx1[count2], yy1[count2], zz1[count2], xx1[count6], yy1[count6], zz1[count6]);
              if (disttmp < distmin)
              {
                distmin = disttmp;
                count4  = count3;
              }
            }
          }
          used_tmp_D4R[count4] = 1;
          D4Ratoms[count5][count7][count8] = tmp_D4R[count4];
        }
      }

    }

    // checking  D4Ratoms
    for (count1 = 0; count1 <= 6; count1++)
    {
      printf(" ms4b D4R %d\n", count1);
      for (count2 = 0; count2 <= 8; count2++)
      {
        printf(" ms4c D4R %d vertex and atoms ", count2);
        for (count3 = 0; count3 <= 8; count3++)
        {
          printf(" %d ",D4Ratoms[count1][count2][count3]);
        }
        printf("\n");
      }
    }





    InpF = fopen(main_files_nameF,"r");
    for (count3 = 1; count3 <= number_input_files; count3++)
    {
      fgets(line, maxline, InpF);
      strcpy(s3 , grepnword(line, 0, 1));
      sscanf(s3,"%lf", &energy_config[count3]);
      strcpy(s3 , grepnword(line, 0, 3+number_Ge));
      sscanf(s3,"%s", heteroatom_input_FName[count3]);
      printf(" ms5a %.8f %s\n", energy_config[count3], heteroatom_input_FName[count3]);
      //sscanf(line,"%s %d %s", heteroatom_input_FName[count3], &CoorFType_heteroatom_input[count3], main_Output_FName[count3]); 
    }
    close(InpF);

    OutpF = fopen(main_Output_FName, "w");
    // doing the MainJob
    for (MainCount = 1; MainCount <= number_input_files; MainCount++)
    {
      //strcpy(s1, heteroatom_input_FName[MainCount]);
      //strcpy(s2, "stw_scan_xx_xxxx/xxxxx/conp/");
      //for (count5 = 9; count5 <= 15; count5++)
      //{
      //  s2[count5] = s1[count5];
      //}
      //for (count5 = 17; count5 <= 21; count5++)
      //{
      //  s2[count5] = s1[count5];
      //}
      strcpy(s2, heteroatom_input_FName[MainCount]);
      strcat(s2, ".cif");
      printf(" ms5 %.8f %s %s\n", energy_config[MainCount], heteroatom_input_FName[MainCount], s2);
      strcpy(currentFname, s2);
// aaa
      count5 = 5;
      count4 = getatomcoor(s2, count5);

      Tatomnumber = (atomnumber - 6) / 3;
      printf(" ms5a %d\n", count4);
      for (count5 = 0; count5 < Tatomnumber; count5++)
      {
        disttmp = atomsdist(xatom[count5], yatom[count5], zatom[count5], xx1[count5], yy1[count5], zz1[count5]);
        if (disttmp > 1.0)
        {
          printf(" ms6 SiGe atom ordering mistmacht at number %2d in file %s\n", count5, heteroatom_input_FName[MainCount]);
        }
      }

      // initializing Ge_in_D4R
      for (count1 = 0; count1 < N; count1++)
      {
        for (count2 = 0; count2 < N; count2++)
        {
          Ge_in_D4R[count1][count2] = 0;
        }
      }

      count6 = 0;
      for (count5 = 0; count5 < Tatomnumber; count5++)
      {
        strcpy(s1, atomname[count5]);
        if ( ((s1[0] == 'G') || (s1[0] == 'g')) && ((s1[1] == 'E') || (s1[1] == 'e')) )
        {
          count6++;
          Ge_atom[count6][0] = count5;     // OJO [ count6 ] Ge index, [ =0 itself, =1 D4R (1-6) & =2 vertex position site inside D4R (1-8) ]
          if (count5>47)
          {
            Ge_atom[count6][1] = 7;
          }
          for (count1 = 1; count1 <= 6; count1++)
          {
            for (count2 = 1; count2 <= 8; count2++)
            {
              if ( D4Ratoms[count1][count2][1] == count5)
              {
                Ge_atom[count6][1] = count1;
                Ge_atom[count6][2] = count2;
              }
            }

          }
        }
      } 
      if (count6 != number_Ge)
      {
        printf(" ms7 bad reading Ge atoms in file %s\n", heteroatom_input_FName[MainCount]);
      } 

      printf(" ms8 Ge in file %s at sites and in in D4R and vertex\n", heteroatom_input_FName[MainCount]);
      for (count5 = 1; count5 <= number_Ge; count5++)
      {
        printf("running_Ge %2d Ge_atom %2d in D4R %d at vertext %d\n", count5, Ge_atom[count5][0], Ge_atom[count5][1], Ge_atom[count5][2]);  
      }
      printf("\n");

      // processing Ge at each D4R
      for (count1 = 1; count1 <= number_Ge; count1++)
      {
        count2 = Ge_atom[count1][1];
        count7 = Ge_in_D4R[count2][0];
        count7++;
        Ge_in_D4R[count2][0] = count7;
        count4 = Ge_atom[count1][0];
        Ge_in_D4R[count2][count7] = count4;
      }
      
      for (count1 = 1; count1 <= 6; count1++)
      {
        count2 = Ge_in_D4R[count1][0];
        printf(" ms9a in D4R %d number of Ge %d\n", count1, count2);
        printf(" ms9b atoms: ");
        for (count7 = 1; count7 <= count2; count7++)
        {
          printf("%d ", Ge_in_D4R[count1][count7]);
        }
        printf("\n");
      }
      printf("\n"); 

      // processing final output
      for (count1 = 1; count1 < 30; count1++)
      {
        D4R_stat[count1] = 0;
      }
      // D4R_stat 
      // 0  number of Ge 0 by D4R
      // 1  number of Ge 1 by D4R
      // 2  number of Ge 2 by D4R
      // 3  number of Ge 2 by D4R as direct neighbour 
      // 4  number of Ge 2 by D4R as face diagonal neighbour
      // 5  number of Ge 2 by D4R as cube diagonal neighbour
      // 6  number of Ge 3 by D4R
      // 7  number of Ge 3 by D4R as 3 atoms chain
      // 8  number of Ge 3 by D4R as 2 direct neighbour and 1 not direct neighbour
      // 9  number of Ge 3 by D4R as 3 face diagonal neighbour
      // 10 number of Ge 4 by D4R 
      // 11 number of Ge 4 by D4R as 1 central atom and 3 direct neighbours
      // 12 number of Ge 4 by D4R as forming a square
      // 13 number of Ge 4 by D4R as 4 atoms chain
      // 14 number of Ge 4 by D4R as 3 atoms chain and 1 isolated
      // 15 number of Ge 4 by D4R as 2 pairs of 2 direct neighbour
      // 16 number of Ge 5 by D4R
      // 17 number of Ge 6 by D4R
      // 18 number of Ge 7 by D4R
      // 19 number of Ge 8 by D4R
      // 20 number of overall isolated Ge atoms
      // 21 numbre of isolated Ge atoms in D4R
      // 22 number of isolated Ge atoms outside D4R, i.e. at T5

      // 23 number of Ge 5 by D4R as forming a square and 1 isolated
      // 24 number of Ge 5 by D4R as 1 central atom and 3 direct neighbours and 1 isolated
      // 25 number of Ge 5 by D4R as 5 atoms chain

      // 26 number of Ge 6 by D4R as forming a square and 2 direct neighbours
      // 27 number of Ge 6 by D4R as forming a square and 2 Ge in face diagonal
      // 28 number of Ge 6 by D4R as 6 atoms chain

      // analyzing pure silica D4R
      D4R_stat[0] = 6;
      tmp_D4R[1]  = 0; tmp_D4R[2]  = 0; tmp_D4R[3]  = 0; tmp_D4R[4]  = 0; tmp_D4R[5]  = 0; tmp_D4R[6]  = 0; tmp_D4R[7]  = 0;
      for (count1 = 1; count1 <= number_Ge; count1++)
      {
        count7 = Ge_atom[count1][1];
        count8 = tmp_D4R[count7];
        count8++;
        tmp_D4R[count7] = count8;
      }
      for (count1 = 1; count1 <= 6; count1++)
      {
        if (tmp_D4R[count1] > 0)
        {
          count8 = D4R_stat[0];
          count8--;
          D4R_stat[0] = count8;
        }
      }  

      // analyzing isolated Ge atoms
      tmp_D4R[1]  = number_Ge;
      tmp_D4R[3]  = 0;
      for (count1 = 1; count1 <= number_Ge; count1++)
      {
        if (Ge_atom[count1][0] > 47)
        {
          count8 = tmp_D4R[3];
          count8++;
          tmp_D4R[3] = count8;
        }
      }
      count8 = tmp_D4R[3];
      tmp_D4R[2]  = number_Ge - count8;
      // overall isolated Ge atoms, isolated Ge atoms in D4R, isolated Ge atoms outside D4R
      for (count1 = 1; count1 <= number_Ge; count1++)
      {
        count3 = Ge_atom[count1][0];
        for (count2 = 1; count2 <= number_Ge; count2++)
        {
          count4 = Ge_atom[count2][0];
          tmp_dist[1] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count4], yy1[count4], zz1[count4]);
          if ( (tmp_dist[1] < 3.8) && (tmp_dist[1] > 0.2) )
          {
            count8 = tmp_D4R[1];
            count8--;
            tmp_D4R[1] = count8;
            if (Ge_atom[count1][1] < 7)
            {
              count8 = tmp_D4R[2];
              count8--;
              tmp_D4R[2] = count8;
            }
            else
            {
              count8 = tmp_D4R[3];
              count8--;
              tmp_D4R[3] = count8;
            }
            count2 = number_Ge + 1;
          }
        }
      }
      D4R_stat[20] = tmp_D4R[1]; 
      D4R_stat[21] = tmp_D4R[2];
      D4R_stat[22] = tmp_D4R[3];

      for (count1 = 1; count1 <= 6; count1++)
      {
        count2 = Ge_in_D4R[count1][0];
        switch (count2)
        {
          case 1:
            count7 = D4R_stat[1];
            count7++;
            D4R_stat[1] = count7;
          break;

          case 2:
            count7 = D4R_stat[2];
            count7++;
            D4R_stat[2] = count7;
            count3 = Ge_in_D4R[count1][1];
            count4 = Ge_in_D4R[count1][2]; 
            //printf(" ms10b %d %d 2 Ge in the %d D4R\n", count3, count4, count1); 
            for (MainCount1 = 1; MainCount1 <= 6; MainCount1++)
            {
              for (MainCount2 = 1; MainCount2 <= 8; MainCount2++)
              {  
                //count8 = D4Ratoms[MainCount1][MainCount2][1];
                //printf(" ms10c %d running atom compared to %d at %d D4R %d column\n", count8, count3, MainCount1, MainCount2);
                //if (count8 == count3)
                if (D4Ratoms[MainCount1][MainCount2][1] == count3)
                {
                  //printf(" ms10d %d 1st Ge in the %d D4R %d column\n", count3, count1, MainCount2);
                  for (MainCount3 = 2; MainCount3 <= 8; MainCount3++)
                  {
                    if (D4Ratoms[MainCount1][MainCount2][MainCount3] == count4)
                    {
                      //printf(" ms10e %d %d 2 Ge in the %d D4R %d column %d row\n", count3, count4, count1, MainCount2, MainCount3);
                      if ((MainCount3>1) && (MainCount3<5))
                      {
                        count7 = D4R_stat[3];
                        count7++;
                        D4R_stat[3] = count7;
                      }
                      else
                      {
                        if ((MainCount3>4) && (MainCount3<8))
                        {
                           count7 = D4R_stat[4];
                           count7++;
                           D4R_stat[4] = count7;
                        }
                        else
                        {
                          count7 = D4R_stat[5];
                          count7++;
                          D4R_stat[5] = count7;
                        }
                      }
                      MainCount3 = 9;
                    }
                    //MainCount3 = 9;
                  }
                  MainCount2 = 9; MainCount1 = 7;
                }
                //MainCount2 = 9; MainCount1 = 7;
              }
            }
          break;

          case 3:
            count7 = D4R_stat[6];
            count7++;
            D4R_stat[6] = count7; 

            count3 = Ge_in_D4R[count1][1];
            count4 = Ge_in_D4R[count1][2]; 
            count5 = Ge_in_D4R[count1][3];
            tmp_dist[1] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count4], yy1[count4], zz1[count4]);
            tmp_dist[2] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count5], yy1[count5], zz1[count5]);
            tmp_dist[3] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count5], yy1[count5], zz1[count5]);            
            tmp_D4R[1]  = 0; tmp_D4R[2]  = 0; tmp_D4R[3]  = 0;
            for (MainCount1 = 1; MainCount1 <= 3; MainCount1++)
            {
              if (tmp_dist[MainCount1] < 4.0) 
              {
                count8 = tmp_D4R[1];
                count8++;
                tmp_D4R[1] = count8;
              }
              else
              {
                if (tmp_dist[MainCount1] < 5.0)
                {
                  count8 = tmp_D4R[2];
                  count8++;
                  tmp_D4R[2] = count8;
                }
                else
                {
                  count8 = tmp_D4R[3];
                  count8++;
                  tmp_D4R[3] = count8;
                }
              }
            }
            if ( (tmp_D4R[1] == 2) && (tmp_D4R[2] == 1) )
            {
              count7 = D4R_stat[7];
              count7++;
              D4R_stat[7] = count7;
            } 
            else
            {
              if ( tmp_D4R[2] == 3 )
              {
                count7 = D4R_stat[9];
                count7++;
                D4R_stat[9] = count7;
              }
              else
              {
                count7 = D4R_stat[8];
                count7++;
                D4R_stat[8] = count7;
              }
            }

          break;

          case 4:
            count7 = D4R_stat[10];
            count7++;
            D4R_stat[10] = count7;

            count3 = Ge_in_D4R[count1][1];
            count4 = Ge_in_D4R[count1][2]; 
            count5 = Ge_in_D4R[count1][3];
            count6 = Ge_in_D4R[count1][4];
            tmp_dist[1] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count4], yy1[count4], zz1[count4]);
            tmp_dist[2] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count5], yy1[count5], zz1[count5]);
            tmp_dist[3] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[4] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count5], yy1[count5], zz1[count5]);            
            tmp_dist[5] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[6] = atomsdist(xx1[count5], yy1[count5], zz1[count5], xx1[count6], yy1[count6], zz1[count6]);
            tmp_D4R[1]  = 0; tmp_D4R[2]  = 0; tmp_D4R[3]  = 0;
            for (MainCount1 = 1; MainCount1 <= 6; MainCount1++)
            {
              if (tmp_dist[MainCount1] < 4.0)              // neighbour
              {
                count8 = tmp_D4R[1];
                count8++;
                tmp_D4R[1] = count8;
              }
              else
              {
                if (tmp_dist[MainCount1] < 5.0)            // face diagonal
                {
                  count8 = tmp_D4R[2];
                  count8++;
                  tmp_D4R[2] = count8;
                }
                else
                {                                          // cube diagonal
                  count8 = tmp_D4R[3];
                  count8++;
                  tmp_D4R[3] = count8;
                }
              }
            }
            if ( (tmp_D4R[1] == 3) && (tmp_D4R[2] == 3) )
            {
              count7 = D4R_stat[11];
              count7++;
              D4R_stat[11] = count7;
            } 
            if ( (tmp_D4R[1] == 4) && (tmp_D4R[2] == 2) )
            {
              count7 = D4R_stat[12];
              count7++;
              D4R_stat[12] = count7;
            }
            if ( (tmp_D4R[1] == 3) && (tmp_D4R[2] == 2) && (tmp_D4R[3] == 1))
            {
              count7 = D4R_stat[13];
              count7++;
              D4R_stat[13] = count7;
            }
            if ( (tmp_D4R[1] == 2) && (tmp_D4R[2] == 3) && (tmp_D4R[3] == 1))
            {
              count7 = D4R_stat[14];
              count7++;
              D4R_stat[14] = count7;
            }
            if ( (tmp_D4R[1] == 2) && (tmp_D4R[2] == 2) && (tmp_D4R[3] == 2))
            {
              count7 = D4R_stat[15];
              count7++;
              D4R_stat[15] = count7;
            }
          break;
        
          case 5: 
            count7 = D4R_stat[16];
            count7++;
            D4R_stat[16] = count7;

            count3 = Ge_in_D4R[count1][1];
            count4 = Ge_in_D4R[count1][2]; 
            count5 = Ge_in_D4R[count1][3];
            count6 = Ge_in_D4R[count1][4];
            count7 = Ge_in_D4R[count1][5];
            tmp_dist[1] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count4], yy1[count4], zz1[count4]);
            tmp_dist[2] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count5], yy1[count5], zz1[count5]);
            tmp_dist[3] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[4] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count7], yy1[count7], zz1[count7]);
            tmp_dist[5] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count5], yy1[count5], zz1[count5]);            
            tmp_dist[6] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[7] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count7], yy1[count7], zz1[count7]);
            tmp_dist[8] = atomsdist(xx1[count5], yy1[count5], zz1[count5], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[9] = atomsdist(xx1[count5], yy1[count5], zz1[count5], xx1[count7], yy1[count7], zz1[count7]);
            tmp_dist[10]= atomsdist(xx1[count6], yy1[count6], zz1[count6], xx1[count7], yy1[count7], zz1[count7]);
            tmp_D4R[1]  = 0; tmp_D4R[2]  = 0; tmp_D4R[3]  = 0;
            for (MainCount1 = 1; MainCount1 <= 10; MainCount1++)
            {
              if (tmp_dist[MainCount1] < 4.0)              // neighbour
              {
                count8 = tmp_D4R[1];
                count8++;
                tmp_D4R[1] = count8;
              }
              else
              {
                if (tmp_dist[MainCount1] < 5.0)            // face diagonal
                {
                  count8 = tmp_D4R[2];
                  count8++;
                  tmp_D4R[2] = count8;
                }
                else
                {                                          // cube diagonal
                  count8 = tmp_D4R[3];
                  count8++;
                  tmp_D4R[3] = count8;
                }
              }
            }
            printf(" D4R %d ms_test_n56 number_neighbour %d number_face_diagonal %d number_cube_diagonal %d\n", count1, tmp_D4R[1], tmp_D4R[2], tmp_D4R[3]);
            if ( (tmp_D4R[1] == 5) && (tmp_D4R[2] == 4) && (tmp_D4R[3] == 1) )
            {
              count7 = D4R_stat[23];
              count7++;
              D4R_stat[23] = count7;
            } 
            if ( (tmp_D4R[1] == 3) && (tmp_D4R[2] == 6) && (tmp_D4R[3] == 1) )
            {
              count7 = D4R_stat[24];
              count7++;
              D4R_stat[24] = count7;
            }
            if ( (tmp_D4R[1] == 4) && (tmp_D4R[2] == 4) && (tmp_D4R[3] == 2) )
            {
              count7 = D4R_stat[25];
              count7++;
              D4R_stat[25] = count7;
            }
          break;

          case 6: 
            count7 = D4R_stat[17];
            count7++;
            D4R_stat[17] = count7;

            count3 = Ge_in_D4R[count1][1];
            count4 = Ge_in_D4R[count1][2]; 
            count5 = Ge_in_D4R[count1][3];
            count6 = Ge_in_D4R[count1][4];
            count7 = Ge_in_D4R[count1][5];
            count8 = Ge_in_D4R[count1][6];
            tmp_dist[1] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count4], yy1[count4], zz1[count4]);
            tmp_dist[2] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count5], yy1[count5], zz1[count5]);
            tmp_dist[3] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[4] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count7], yy1[count7], zz1[count7]);
            tmp_dist[5] = atomsdist(xx1[count3], yy1[count3], zz1[count3], xx1[count8], yy1[count8], zz1[count8]);
            tmp_dist[6] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count5], yy1[count5], zz1[count5]);            
            tmp_dist[7] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[8] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count7], yy1[count7], zz1[count7]);
            tmp_dist[9] = atomsdist(xx1[count4], yy1[count4], zz1[count4], xx1[count8], yy1[count8], zz1[count8]);
            tmp_dist[10]= atomsdist(xx1[count5], yy1[count5], zz1[count5], xx1[count6], yy1[count6], zz1[count6]);
            tmp_dist[11]= atomsdist(xx1[count5], yy1[count5], zz1[count5], xx1[count7], yy1[count7], zz1[count7]);
            tmp_dist[12]= atomsdist(xx1[count5], yy1[count5], zz1[count5], xx1[count8], yy1[count8], zz1[count8]);
            tmp_dist[13]= atomsdist(xx1[count6], yy1[count6], zz1[count6], xx1[count7], yy1[count7], zz1[count7]);
            tmp_dist[14]= atomsdist(xx1[count6], yy1[count6], zz1[count6], xx1[count8], yy1[count8], zz1[count8]);
            tmp_dist[15]= atomsdist(xx1[count7], yy1[count7], zz1[count7], xx1[count8], yy1[count8], zz1[count8]);
            tmp_D4R[1]  = 0; tmp_D4R[2]  = 0; tmp_D4R[3]  = 0;
            for (MainCount1 = 1; MainCount1 <= 15; MainCount1++)
            {
              if (tmp_dist[MainCount1] < 4.0)              // neighbour
              {
                count8 = tmp_D4R[1];
                count8++;
                tmp_D4R[1] = count8;
              }
              else
              {
                if (tmp_dist[MainCount1] < 5.0)            // face diagonal
                {
                  count8 = tmp_D4R[2];
                  count8++;
                  tmp_D4R[2] = count8;
                }
                else
                {                                          // cube diagonal
                  count8 = tmp_D4R[3];
                  count8++;
                  tmp_D4R[3] = count8;
                }
              }
            }
            printf(" D4R %d ms_test_n56 number_neighbour %d number_face_diagonal %d number_cube_diagonal %d\n", count1, tmp_D4R[1], tmp_D4R[2], tmp_D4R[3]);
            if ( (tmp_D4R[1] == 7) && (tmp_D4R[2] == 6) && (tmp_D4R[3] == 2) )
            {
              count7 = D4R_stat[26];
              count7++;
              D4R_stat[26] = count7;
            } 
            if ( (tmp_D4R[1] == 6) && (tmp_D4R[2] == 7) && (tmp_D4R[3] == 2) )
            {
              count7 = D4R_stat[27];
              count7++;
              D4R_stat[27] = count7;
            }
            if ( (tmp_D4R[1] == 6) && (tmp_D4R[2] == 6) && (tmp_D4R[3] == 3) )
            {
              count7 = D4R_stat[28];
              count7++;
              D4R_stat[28] = count7;
            }
          break;

          case 7:
            count7 = D4R_stat[18];
            count7++;
            D4R_stat[18] = count7;
          break;

          case 8:
            count7 = D4R_stat[19];
            count7++;
            D4R_stat[19] = count7;
          break;
        }
      }
      printf(" ms10a %d configuration summary information\n", MainCount);
      printf("%s    // Current File Name\n", currentFname);
      printf("%f   // relative lattice energy\n", energy_config[MainCount]);
      printf("%2d         // 0  number of Ge 0 by D4R\n", D4R_stat[0]);
      printf("%2d         // 1  number of Ge 1 by D4R\n", D4R_stat[1]);
      printf("%2d         // 2  number of Ge 2 by D4R\n", D4R_stat[2]);
      printf("%2d         // 3  number of Ge 2 by D4R as direct neighbour\n", D4R_stat[3]);
      printf("%2d         // 4  number of Ge 2 by D4R as face diagonal neighbour\n", D4R_stat[4]);
      printf("%2d         // 5  number of Ge 2 by D4R as cube diagonal neighbour\n", D4R_stat[5]);
      printf("%2d         // 6  number of Ge 3 by D4R\n", D4R_stat[6]);
      printf("%2d         // 7  number of Ge 3 by D4R as 3 atoms chain\n", D4R_stat[7]);
      printf("%2d         // 8  number of Ge 3 by D4R as 2 direct neighbour and 1 not direct neighbour\n", D4R_stat[8]);
      printf("%2d         // 9  number of Ge 3 by D4R as 3 face diagonal neighbour\n", D4R_stat[9]);
      printf("%2d         // 10 number of Ge 4 by D4R\n", D4R_stat[10]);
      printf("%2d         // 11 number of Ge 4 by D4R as 1 central atom and 3 direct neighbours\n", D4R_stat[11]);
      printf("%2d         // 12 number of Ge 4 by D4R as forming a square\n", D4R_stat[12]);
      printf("%2d         // 13 number of Ge 4 by D4R as 4 atoms chain\n", D4R_stat[13]);
      printf("%2d         // 14 number of Ge 4 by D4R as 3 atoms chain and 1 isolated\n", D4R_stat[14]);
      printf("%2d         // 15 number of Ge 4 by D4R as 2 pairs of 2 direct neighbour\n", D4R_stat[15]);
      printf("%2d         // 16 number of Ge 5 by D4R\n", D4R_stat[16]);
      printf("%2d         // 17 number of Ge 5 by D4R as forming a square and 1 isolated\n", D4R_stat[23]);
      printf("%2d         // 18 number of Ge 5 by D4R as 1 central atom and 3 direct neighbours and 1 isolated\n", D4R_stat[24]);
      printf("%2d         // 19 number of 5 by D4R as 5 atoms chain\n", D4R_stat[25]);
      printf("%2d         // 19 number of Ge 6 by D4R\n", D4R_stat[17]);
      printf("%2d         // 20 number of Ge 6 by D4R as forming a square and 2 direct neighbours\n", D4R_stat[26]);
      printf("%2d         // 20 number of Ge 6 by D4R as forming a square and 2 Ge in face diagonal\n", D4R_stat[27]);
      printf("%2d         // 21 numbre of Ge 6 by D4R as 6 atoms chai\n", D4R_stat[28]);
      printf("%2d         // 22 number of Ge 7 by D4R\n", D4R_stat[18]);
      printf("%2d         // 23 number of Ge 8 by D4R\n", D4R_stat[19]);
      printf("%2d         // 24 number of overall isolated Ge atoms\n", D4R_stat[20]);
      printf("%2d         // 25 numbre of isolated Ge atoms in D4R\n", D4R_stat[21]);
      printf("%2d         // 26 number of isolated Ge atoms outside D4R, i.e. at T5\n", D4R_stat[22]);

      // Ge population by crsytallographic site
      tmp_D4R[1]  = 0; tmp_D4R[2]  = 0; tmp_D4R[3]  = 0; tmp_D4R[4]  = 0; tmp_D4R[5]  = 0;
      //for (count1 = 1; count1 <= number_Ge; count1++)
      for (count1 = 1; count1 <= number_Ge; count1++)
      {
        if (Ge_atom[count1][0] < 12)
        {
          count8 = tmp_D4R[1];
          count8++;
          tmp_D4R[1] = count8;
          printf("ms1x counting T1 site occupancy: running_scan %2d Ge_atom %2d count %2d\n", count1, Ge_atom[count1][0], tmp_D4R[1]);
        }
        if ( (Ge_atom[count1][0] > 11) && (Ge_atom[count1][0] < 24) )
        {
          count8 = tmp_D4R[2];
          count8++;
          tmp_D4R[2] = count8;
        }
        if ( (Ge_atom[count1][0] > 23) && (Ge_atom[count1][0] < 36) )
        {
          count8 = tmp_D4R[3];
          count8++;
          tmp_D4R[3] = count8;
        }
        if ( (Ge_atom[count1][0] > 35) && (Ge_atom[count1][0] < 48) )
        {
          count8 = tmp_D4R[4];
          count8++;
          tmp_D4R[4] = count8;
        }
        if (Ge_atom[count1][0] > 47)
        {
          count8 = tmp_D4R[5];
          count8++;
          tmp_D4R[5] = count8;
        }
      }
      printf("%2d   %.5f %.5f %.5f   // 20 number of Ge at crystallographic site 1 with the given coordenates\n", tmp_D4R[1], xatom[0], yatom[0], zatom[0]);
      printf("%2d   %.5f %.5f %.5f   // 21 number of Ge at crystallographic site 2 with the given coordenates\n", tmp_D4R[2], xatom[12], yatom[12], zatom[12]);
      printf("%2d   %.5f %.5f %.5f   // 22 number of Ge at crystallographic site 3 with the given coordenates\n", tmp_D4R[3], xatom[24], yatom[24], zatom[24]);
      printf("%2d   %.5f %.5f %.5f   // 23 number of Ge at crystallographic site 4 with the given coordenates\n", tmp_D4R[4], xatom[36], yatom[36], zatom[36]);
      printf("%2d   %.5f %.5f %.5f   // 24 number of Ge at crystallographic site 5 with the given coordenates\n", tmp_D4R[5], xatom[48], yatom[48], zatom[48]);
      printf("\n");

      // printing out results to file
      fprintf(OutpF," %d configuration summary information\n", MainCount);
      fprintf(OutpF," %s    // Current File Name\n", currentFname);
      fprintf(OutpF," %f   // relative lattice energy\n", energy_config[MainCount]);
      fprintf(OutpF,"%2d         // 0  number of Ge 0 by D4R\n", D4R_stat[0]);
      fprintf(OutpF,"%2d         // 1  number of Ge 1 by D4R\n", D4R_stat[1]);
      fprintf(OutpF,"%2d         // 2  number of Ge 2 by D4R\n", D4R_stat[2]);
      fprintf(OutpF,"%2d         // 3  number of Ge 2 by D4R as direct neighbour\n", D4R_stat[3]);
      fprintf(OutpF,"%2d         // 4  number of Ge 2 by D4R as face diagonal neighbour\n", D4R_stat[4]);
      fprintf(OutpF,"%2d         // 5  number of Ge 2 by D4R as cube diagonal neighbour\n", D4R_stat[5]);
      fprintf(OutpF,"%2d         // 6  number of Ge 3 by D4R\n", D4R_stat[6]);
      fprintf(OutpF,"%2d         // 7  number of Ge 3 by D4R as 3 atoms chain\n", D4R_stat[7]);
      fprintf(OutpF,"%2d         // 8  number of Ge 3 by D4R as 2 direct neighbour and 1 not direct neighbour\n", D4R_stat[8]);
      fprintf(OutpF,"%2d         // 9  number of Ge 3 by D4R as 3 face diagonal neighbour\n", D4R_stat[9]);
      fprintf(OutpF,"%2d         // 10 number of Ge 4 by D4R\n", D4R_stat[10]);
      fprintf(OutpF,"%2d         // 11 number of Ge 4 by D4R as 1 central atom and 3 direct neighbours\n", D4R_stat[11]);
      fprintf(OutpF,"%2d         // 12 number of Ge 4 by D4R as forming a square\n", D4R_stat[12]);
      fprintf(OutpF,"%2d         // 13 number of Ge 4 by D4R as 4 atoms chain\n", D4R_stat[13]);
      fprintf(OutpF,"%2d         // 14 number of Ge 4 by D4R as 3 atoms chain and 1 isolated\n", D4R_stat[14]);
      fprintf(OutpF,"%2d         // 15 number of Ge 4 by D4R as 2 pairs of 2 direct neighbour\n", D4R_stat[15]);
      fprintf(OutpF,"%2d         // 16 number of Ge 5 by D4R\n", D4R_stat[16]);
      fprintf(OutpF,"%2d         // 17 number of Ge 5 by D4R as forming a square and 1 isolated\n", D4R_stat[23]);
      fprintf(OutpF,"%2d         // 18 number of Ge 5 by D4R as 1 central atom and 3 direct neighbours and 1 isolated\n", D4R_stat[24]);
      fprintf(OutpF,"%2d         // 19 number of 5 by D4R as 5 atoms chain\n", D4R_stat[25]);
      fprintf(OutpF,"%2d         // 19 number of Ge 6 by D4R\n", D4R_stat[17]);
      fprintf(OutpF,"%2d         // 20 number of Ge 6 by D4R as forming a square and 2 direct neighbours\n", D4R_stat[26]);
      fprintf(OutpF,"%2d         // 20 number of Ge 6 by D4R as forming a square and 2 Ge in face diagonal\n", D4R_stat[27]); 
      fprintf(OutpF,"%2d         // 21 numbre of Ge 6 by D4R as 6 atoms chai\n", D4R_stat[28]);
      fprintf(OutpF,"%2d         // 22 number of Ge 7 by D4R\n", D4R_stat[18]);
      fprintf(OutpF,"%2d         // 23 number of Ge 8 by D4R\n", D4R_stat[19]);
      fprintf(OutpF,"%2d         // 24 number of overall isolated Ge atoms\n", D4R_stat[20]);
      fprintf(OutpF,"%2d         // 25 numbre of isolated Ge atoms in D4R\n", D4R_stat[21]);
      fprintf(OutpF,"%2d         // 26 number of isolated Ge atoms outside D4R, i.e. at T5\n", D4R_stat[22]);
      fprintf(OutpF,"%2d   %.5f %.5f %.5f   // 27 number of Ge at crystallographic site 1 with the given coordenates\n", tmp_D4R[1], xatom[0], yatom[0], zatom[0]);
      fprintf(OutpF,"%2d   %.5f %.5f %.5f   // 28 number of Ge at crystallographic site 2 with the given coordenates\n", tmp_D4R[2], xatom[12], yatom[12], zatom[12]);
      fprintf(OutpF,"%2d   %.5f %.5f %.5f   // 29 number of Ge at crystallographic site 3 with the given coordenates\n", tmp_D4R[3], xatom[24], yatom[24], zatom[24]);
      fprintf(OutpF,"%2d   %.5f %.5f %.5f   // 30 number of Ge at crystallographic site 4 with the given coordenates\n", tmp_D4R[4], xatom[36], yatom[36], zatom[36]);
      fprintf(OutpF,"%2d   %.5f %.5f %.5f   // 31 number of Ge at crystallographic site 5 with the given coordenates\n", tmp_D4R[5], xatom[48], yatom[48], zatom[48]);
      fprintf(OutpF,"\n");
//aaa */      
    }

    fclose(OutpF);

}    /***********    end of MainJob     ************/

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
         //for (c1 = 0; c1 <= atomnumber-1; c1++)
         //      {
         //        printf("%s %c %.5f %.5f %.5f\n", atomname[c1], atom_c_s[c1], xatom[c1], yatom[c1], zatom[c1]);  
         //      }

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

x11, y11, z11 are the coordinates of the central atom of the angle
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
/*    printf(" cosin of the angle %f\n", angle1); */
    
    // correcting 180 and 0 degrees for error induced by numerical truncation 
    if ( (angle1<-1.0) &&  (angle1>-1.01) )
    {
      angle1 = -1.0;
    }
    else
    if ( (angle1>1.0) &&  (angle1<1.01) )
    {
      angle1 = 1.0;
    }
	angle1 = acos(angle1)*180.0/M_PI;
/*    printf("              angle %f\n", angle1); */
    
    
    return angle1;
}
/************************** end of atomsangle ********************************/


/*****************************************************************************
this function gives the angles between three atoms in a cartesian non-periodic
system
x11, y11, z11 are the coordinates of the central atom of the angle
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
    
    // correcting 180 and 0 degrees for error induced by numerical truncation 
    if ( (angle1<-1.0) &&  (angle1>-1.01) )
    {
      angle1 = -1.0;
    }
    else
    if ( (angle1>1.0) &&  (angle1<1.01) )
    {
      angle1 = 1.0;
    }
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
        strcpy(grepnword1, ""); memset(grepnword1, 0,  sizeof(grepnword1));
        //printf("in grepnword1 internal_string %s and function_argument %s\n", grepnword1, line1);
        
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

        //printf("in grepnword2 internal_string %s and function_argument %s\n", grepnword1, line1);

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

        return grepnmword2;

} /***********************  end of grepnmword *******************************/
                                                                               

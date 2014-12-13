/*
    agruments: Input pic(garb34.pgm), sig value(1)
    outputs: mag.pgm, peaks.pgm, final.pgm
    terminal prompt for input pct (4.7% is a good value)
*/


#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#include <stdlib.h>
#define  PICSIZE 256
#define  MAXMASK 7 //6sig+1 

         int    pic[PICSIZE][PICSIZE];
         double outpicx[PICSIZE][PICSIZE];
         double outpicy[PICSIZE][PICSIZE];
         double maskx[MAXMASK][MAXMASK]; 
         double masky[MAXMASK][MAXMASK];
         double mag[PICSIZE][PICSIZE];
         double peaks[PICSIZE][PICSIZE];
         double final[PICSIZE][PICSIZE];
         int freq[PICSIZE]; // for auto thresholds 

main(argc,argv)
int argc;
char **argv;
{
        int     i,j,x,y,mr,centx,centy;
        double  maskval,sum,sum1,sum2,sig,maxival,slope,cutoff,pct;
        double hi, lo; //for double thresholding, might need ints
        FILE    *fo1, *fo2, *fo3,*fp1, *fopen();
        char    *foobar;
        char throwaway[80]; // throw away the pgm header info
        int MTD;

        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");
     
        
        argc--; argv++;
        foobar = *argv;
        sig = atof(foobar);
        
        printf("Enter an input percentage: \n");
        scanf("%lf",&pct);    

        fo1=fopen("mag.pgm","wb"); // open up the outputfile
        fo2=fopen("peaks.pgm","wb"); // peaks output file
        fo3=fopen("final.pgm","wb"); // peaks output file


        // throw away the 3 initial formatting lines
        fgets(throwaway,80,fp1);
        fgets(throwaway,80,fp1);
        fgets(throwaway,80,fp1);

        mr = (int)(sig * 3);
        centx = (MAXMASK / 2);
        centy = (MAXMASK / 2);

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                }
        }
        // changed to first X-derivative 
        for (x=-mr;x<=mr;x++)
        {  for (y=-mr;y<=mr;y++)
           {
              maskval = ((-1*(x/(2*M_PI*sig*sig*sig*sig)))*
                      (exp(-1*(((x*x)+(y*y))/(2*(sig*sig))))));              
              (masky[x+centy][y+centx]) = maskval; // correct
           }
        }

        //adding Y-derivative
        for (x=-mr;x<=mr;x++)
        {  for (y=-mr;y<=mr;y++)
           {
              maskval = ((-1*(y/(2*M_PI*sig*sig*sig*sig)))*
                      (exp(-1*(((x*x)+(y*y))/(2*(sig*sig))))));              
              (maskx[x+centy][y+centx]) = maskval; // correct
           }
        }
        // convolution, changed to double-up
        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sum1 = 0;
             sum2 = 0;
             for (x=-mr;x<=mr;x++)
             {
                for (y=-mr;y<=mr;y++)
                {
                   sum1 += pic[i+x][j+y] * maskx[x+centy][y+centx];
                   sum2 += pic[i+x][j+y] * masky[x+centy][y+centx];
                }
             }
             outpicx[i][j] = sum1; // conv x
             outpicy[i][j] = sum2; //conv y
             //conv[i][j] = sum; // this was here before
          }
        }

       // bring in the sqrt code from sobel below here

        maxival = 0;
        for (i=mr;i<256-mr;i++)// magnitude formula
        { for (j=mr;j<256-mr;j++)
          {
             mag[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                      (outpicy[i][j]*outpicy[i][j])));
             if (mag[i][j] > maxival)
                maxival = mag[i][j];

           }
        }

         //pgm formatting        
        fprintf(fo1,"P5\n");
        fprintf(fo1,"256 256\n");
        fprintf(fo1,"255\n");

        // just get the scaled values
        // print to first output image
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             mag[i][j] = (mag[i][j] / maxival) * 255;            
             fprintf(fo1,"%c",(char)((int)(mag[i][j])));
             
            }
          }



        // peaks code from slide
      for(i=mr;i<256-mr;i++)
      {
        for(j=mr;j<256-mr;j++)
        {    

          if((outpicx[i][j]) == 0.0) 
          {
            outpicx[i][j] = .00001;
          }

          slope = outpicy[i][j]/outpicx[i][j];
          if( (slope <= .4142)&&(slope > -.4142))
          {
            if((mag[i][j] > mag[i][j-1])&&(mag[i][j] > mag[i][j+1]))
            {
            peaks[i][j] = 255;
            }
          }
          else if( (slope <= 2.4142)&&(slope > .4142))
          {
            if((mag[i][j] > mag[i-1][j-1])&&(mag[i][j] > mag[i+1][j+1]))
            {
             peaks[i][j] = 255;
            }
           }
          else if( (slope <= -.4142)&&(slope > -2.4142))
          {
            if((mag[i][j] > mag[i+1][j-1])&&(mag[i][j] > mag[i-1][j+1]))
            {
             peaks[i][j] = 255;
            }
          }
          else
          {
            if((mag[i][j] > mag[i-1][j])&&(mag[i][j] > mag[i+1][j]))
            {
             peaks[i][j] = 255;
            }
          }
        }
      }


        //pgm formatting        
        fprintf(fo2,"P5\n");
        fprintf(fo2,"256 256\n");
        fprintf(fo2,"255\n");

        //output peaks to file
        for (i=0;i<256;i++)
        {
            for (j=0;j<256;j++)
            {
                         
             fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
             
            }
         }

         // begin selecting thresholds

           //initiliaze frequencies to 0
        for(i=0;i<PICSIZE;i++)
        {          
          freq[i]=0;          
        }
          //build up the histogram
        for(i=0;i<PICSIZE;i++)
        {
          for(j=0;j<PICSIZE;j++)
          {              
            freq[(int)(mag[i][j])]++;
          }

        }

        cutoff=pct*PICSIZE*PICSIZE;
        sum=0;
        for(j=255;j>0;j--)
        {
          sum += freq[j];
          if(sum>cutoff)
            {
               hi=j;
               break;
            }
        }

        lo=.35*hi;


        // begin double thresholding 

        // lo= 30;
        // hi= 160;

       
        for(i=0;i<PICSIZE;i++)
        {
          for(j=0;j<PICSIZE;j++)
          {
            if(peaks[i][j]==255)
            {
              if(mag[i][j]>hi)
              {
                peaks[i][j]=0; //off
                final[i][j]= 255; //on
              }
              else if (mag[i][j]<lo) 
              {
                peaks[i][j]=0; 
                final[i][j]=0;
                             
              }
              else
              {
                // pall the mid stuff
              }              
            }
          
          }// end outter j for
        }//end outter i for

        MTD =1; // more to do
              while(MTD==1)
              {
                MTD=0;
                for(i=0;i<PICSIZE;i++)
                {
                  for(j=0;j<PICSIZE;j++)
                  {
                    if(peaks[i][j]==255)
                    {
                      for(x=-1;x<2;x++)
                      {
                        for(y=-1;y<2;y++)
                        {
                          if(final[i+x][j+y]==255)
                          {
                            peaks[i][j]=0;
                            final[i][j]=255;
                            MTD=1; // still more to do
                          }
                          
                          
                        }
                      }
                    }
                     
                  }
                }

              } 

        //pgm formatting        
        fprintf(fo3,"P5\n");
        fprintf(fo3,"256 256\n");
        fprintf(fo3,"255\n");

        //output peaks to file
        for (i=0;i<256;i++)
        {
            for (j=0;j<256;j++)
            {
                         
             fprintf(fo3,"%c",(char)((int)(final[i][j])));
             
            }
         }
 

}//end main



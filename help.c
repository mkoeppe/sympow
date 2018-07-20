#include "sympow.h"

void help_message()
{printf("\nsympow %s takes options:\n",VERSION);
 printf(" -bound #      an upper BOUND for how many ap to compute\n");
 printf(" -help         print the help message and exit\n");
 printf(" -info [] []   only report local information for primes/sympows\n");
 printf("               1st argument is prime range, 2nd is sympow range\n");
 printf(" -local        only report local information (bad primes)\n");
 printf(" -curve []     input a curve in [a1,a2,a3,a4,a6] form\n");
 printf(" -label []     get a label to the given curve\n");
 printf(" -quiet        turn off some messages\n");
 printf(" -verbose      turn on some messages\n");
 printf(" -rootno #     compute the root number of the #th symmetric power\n");
 printf(" -moddeg       compute the modular degree\n");
 printf(" -analrank     compute the analytic rank\n");
 printf(" -sloppy []    for use with -analrank; have X sloppy digits\n");
 printf(" -nocm         abort if curve has complex multiplication\n");
 printf(" -noqt         ignore even powers of non-minimal quad twists\n");
 printf(" -noqdcheck    don't check if quad-double stuff works\n");
 printf(" -mdspeed []   speed for moddeg; 2.0 is default, 0.0 is proof\n");
 printf(" -hecke        compute Hecke symmetric powers for a CM curve\n");
 printf(" -maxtable     set the max size of factor tables: 2^27 default\n" );
 printf(" -sp  []       argument to specify which powers\n");
 printf("               this is a comma separated list\n");
 printf("               in each entry, the 1st datum is the sympow\n");
 printf("               then could come b which turns Bloch-Kato on\n");
 printf("               then could come w# which specifies how many tests\n");
 printf("               then could come s# which says # sloppy digits\n");
 printf("               then must come p# which specifices the precision\n");
 printf("                    or P# which says ignore BOUND for this power\n");
 printf("               then must come d# which says the derivative bound\n");
 printf("                    or D# which says do only this derivative\n");
 printf("                    (neither need be indicated for even powers)\n");
 printf("               default is 2w3s1p32,3bp16d1,4p8\n");
 printf(" -new_data []  will compute inverse Mellin transform mesh for\n");
 printf("               the given data: the format is [sp]d[dv]{h,c}\n");
 printf("               sp is the symmetric power, dv is the derivative,\n");
 printf("               h indicates Hecke powers, and c indicates CM case\n");
 printf("               d[dv] is given only for odd or Hecke powers\n");
 printf("               Examples: 1d3 2 2d1h 3d2 4 4c 5d0 6 7d0h 11d1 12c\n");
 printf("               NOTE: new_data runs a shell script that uses GP\n");
 printf("Other options are used internally/recursively by -new_data\n\n");
 exit(0);}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int create_main(int argc, char **argv);
int view_main(int argc, char **argv);
int denoise_main(int argc, char **argv);
int histo_main(int argc, char **argv);
int htest_main(int argc, char **argv);
int dmodel_main(int argc, char **argv);
int est_dmodel_main(int argc, char **argv);
int lrdmodel_main(int argc, char **argv);

static int usage(char **argv)
{
  printf("\nUsage: %s <command> [options]\n\n", argv[0]);
  printf("Commands:\n");
  printf("       create      -   Create bindump from bam file that subsequent\n");
  printf("                       analysis will take as input.\n");
  printf("       view        -   Inspect contents of bindump and apply filters\n");
  printf("                       to it.\n");
  printf("       histo       -   Show ASCII histogram of frequencies in bindump.\n");
  printf("                       \n");
  printf("       denoise     -   Use the GMMU model to detect and remove a uniform\n");
  printf("                       baseline from the histogram.\n");
  printf("       histotest   -   Use a simple, linear regression based test against\n");
  printf("                       the three fixed models.\n");
  printf("       modeltest   -   Use the GMM based test against the three fixed\n");
  printf("                       models.\n");
  printf("       estmodel    -   Use the GMM to estimate the free model.\n");
  printf("                       \n");
  printf("       lrdmodel    -   Use the GMM to combine the fixed and free models\n");
  printf("                       and assess the delta-log-likelihood.\n\n");

  printf("Please run any of these commands without arguments for usage instructions.\n\n");

  printf("For further documentation please also refer to the online manual\n");
  printf("at https://github.com/clwgg/nQuire or the README file, as well as\nthe publication: http://biorxiv.org/content/early/2017/05/29/143537\n\n");
  return 1;
}

int main(int argc, char **argv)
{
  if (argc < 2) return usage(argv);

  int ret = 0;

  if (strcmp(argv[1], "create") == 0) ret = create_main(argc-1, argv+1);
  else if (strcmp(argv[1], "view") == 0) ret = view_main(argc-1, argv+1);
  else if (strcmp(argv[1], "denoise") == 0) ret = denoise_main(argc-1, argv+1);
  else if (strcmp(argv[1], "histo") == 0) ret = histo_main(argc-1, argv+1);
  else if (strcmp(argv[1], "histotest") == 0) ret = htest_main(argc-1, argv+1);
  else if (strcmp(argv[1], "modeltest") == 0) ret = dmodel_main(argc-1, argv+1);
  else if (strcmp(argv[1], "estmodel") == 0) ret = est_dmodel_main(argc-1, argv+1);
  else if (strcmp(argv[1], "lrdmodel") == 0) ret = lrdmodel_main(argc-1, argv+1);
  else {
    fprintf(stderr, "Unknown command %s\n", argv[1]);
    usage(argv);
  }

  return ret;
}

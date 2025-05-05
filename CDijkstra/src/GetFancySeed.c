#include <sys/types.h>
#include <unistd.h>

unsigned int GetFancySeed(int trulyRandom)
{
    unsigned int seed = 0;
    const char *cmd = "(hostname -i || hostname | sum | awk '{srand(); printf \"%d.%d.%d.%d\n\",256*rand(),256*rand(),$1/256,$1%256}') 2>/dev/null | awk '{for(i=1;i<=NF;i++)if(match($i,\"^[0-9]*\\\\.[0-9]*\\\\.[0-9]*\\\\.[0-9]*$\")){IP=$i;exit}}END{if(!IP)IP=\"127.0.0.1\"; print IP}'";
    FILE *fp=popen(cmd,"r");
    int i, ip[4], host_ip=0;
    if(4!=fscanf(fp," %d.%d.%d.%d ", ip, ip+1, ip+2, ip+3)) fprintf(stderr, "Attempt to get IPv4 address failed:\n%s\n",cmd), exit(1);
    pclose(fp);
    for(i=0;i<4;i++) host_ip = 256*host_ip + ip[i];
    unsigned int dev_random=0;
    if(trulyRandom) {
	fp = fopen("/dev/random","r");
	if(fp){
	    if(1 != fread(&dev_random, sizeof(dev_random),1, fp)) fprintf(stderr, "dev/random read failed\n"), exit(1);
	    fclose(fp);
	}
	else dev_random = lrand48(); // cheap substitute
    }
    //extern int getpid(void), getppid(void);
    seed = host_ip + time(0) + getppid() + getpid() + dev_random;
#if 0
    fprintf(stderr,"%s\n",cmd);
    fprintf(stderr,"%d.%d.%d.%d\n",ip[0],ip[1],ip[2],ip[3]);
    fprintf(stderr,"seed is %ud\n",seed);
#endif
    return seed;
}

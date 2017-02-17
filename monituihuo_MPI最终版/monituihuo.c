#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <sys/time.h>
#include <mpi.h>


#define ZUIDADIANSHU 100	       //最大的点数
#define QISHIDIAN 4	              //起始点
#define CHUSHIWENDU 1400000	      //初始温度
#define WENDUXISHU 0.99	             //降温系数
#define T_MIN 0.000001              //需要降到的最低温

//全局变量
char point_list[ZUIDADIANSHU][20];			//存储点名称
double point_coordinate[ZUIDADIANSHU][2];	      //存储点x y坐标

double P(int *i,int *j,double t,int point_number);        //比较两个路径序列，若长度i>j,返回1.0，否则使用m准则
double distance_sum(int *x,int point_number);             //放回路径序列x的路径长度
void Neighbour(int *father,int *result,int point_number); //产生新的序列
double distance(int x,int y);                             //返回两个点的距离长度
double random0_1(void);                                   //产生随机数,范围在0-1之间
int *random2(int point_number);                           //产生一个以0开头的 0-point_number的不重复序列
double mytime();                                          //返z回当前时间

int main()
{

	int my_rank,comm_sz;              //获取进程号
	double tstart,tstop;              //开始和结束时间
	char filename[81];                //文件名
	double x_position,y_position;    //点的x，y坐标
	FILE *pfile;                    //文件
	int *sequence;                //初始序列
	int point_number;	     //点数量  可以考虑用宏定义定义最大数量
	int i2;			    //for循环专用变量
	int *endlist;              //最短路径序列
	double sumlen;             //最短路径长度
	int **buf2;               //收集各个进程的局部最优解
	int *buf ;              //收集各个进程的局部最优解
	int	*temp ;         //用于交换
	double random ;         //两个路径序列的比较结果
	int *local_i ;          //当前最优解
	int *local_j ;          //新的路径序列
	int c,m,n,k;
	double local_len,nextlen,len,d_len;
 	double t , rate; //初始温度，降温系数
   		              

  MPI_Init(NULL,NULL);               //初始化
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); //获取进程号
  MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);

	tstart=mytime();
	strcpy(filename,"100.txt");

    if((pfile = fopen(filename,"r"))==NULL)
    {
	  printf("Not find %s!\n",filename);
	  return 0;
    }

     /*读取目标文件城市信息*/
    fscanf(pfile,"%d",&point_number);
    if(my_rank==0){
	  printf("numbers 0f city are=%d\n",point_number);
    }

	 for(i2=0;i2<point_number;i2++)
	{
		fscanf(pfile,"%s\t%lf\t%lf",&point_list[i2],&x_position,&y_position);
		point_coordinate[i2][0] = x_position;
		point_coordinate[i2][1] = y_position;
	}

	fclose(pfile);
	srand( (unsigned)time( NULL ) );    //设置种子 初始化随机发生器
	sequence = random2(point_number);

  // MPI_Barrier(MPI_COMM_WORLD);

  /*一些变量赋值 */
	random = 0;
	local_i = malloc(point_number*sizeof(int));
	local_j = malloc(point_number*sizeof(int));
    nextlen=0.0;len=0.0;d_len;
 	t = CHUSHIWENDU;rate = WENDUXISHU;

	for(m=0;m<point_number;m++)
	{
		local_i[m] = sequence[m];
		local_j[m] = 0;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	while(t>T_MIN)
	{
	  do
	  {
		for(c=0;c<20000/comm_sz;c++)
		{
		  Neighbour(local_i,local_j,point_number);
				random = P(local_i,local_j,t,point_number);
				if ( (random == 1.0) || (random > random0_1()) )
				 {
					temp = local_i;
					local_i = local_j;
					local_j = temp;  //保证local_i序列始终是当前起始序列
				 }
		}

		local_len=distance_sum(local_i,point_number);
		if(my_rank==0) nextlen=len;
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Reduce(&local_len,&len,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);  //聚集起来，选取最小的付给len

		if(my_rank==0)
			d_len=len-nextlen;

		MPI_Bcast(&d_len,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	   }while(d_len>0.1);

	   t*=rate;
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	 /* 一些变量赋值*/
	if(my_rank==0)
	{
	   buf2 = malloc(comm_sz*sizeof(int*));
	   for(n=0;n<comm_sz;n++)
	   {
		buf2[n]=(int *)malloc(point_number*sizeof(int));
	   }
	   buf = (int *)malloc(comm_sz*point_number*sizeof(int));

	}

	 //MPI_Barrier(MPI_COMM_WORLD);
	 MPI_Gather(local_i,point_number,MPI_INT,buf,point_number,MPI_INT,0,MPI_COMM_WORLD);    //把各个进程的局部最优解收集起来

	  /*在各个局部最优解中选出最优解*/
	 if(my_rank==0)
	 {
	   for(n=0;n<comm_sz;n++)
	   {
			for(k=0;k<point_number;k++)
			{
			  buf2[n][k]=buf[n*point_number+k];
			}
	   }

	   for(n=1;n<comm_sz;n++)
	   {
		 random = P(buf2[0],buf2[n],t,point_number);
			 if ( (random == 1.0) )
			 {
				temp = buf2[0];
				buf2[0] = buf2[n];
				buf2[n] = temp;    //保证buf2[0]序列始终是当前起始序列
			 }
	   }

	   endlist=buf2[0];
	   printf("last: ");

	   for(i2=0;i2<point_number;i2++)
	   {
			printf("%d ",endlist[i2]);
	   }

	   sumlen = distance_sum(endlist,point_number);
	   printf("\n sumlen :%lf\n",sumlen);
	   free(endlist);
	   free(sequence);
	   tstop=mytime();
	   printf("spend time=%lf\n",tstop-tstart);
	 }

	  MPI_Finalize();
	  return 0;
	}


//产生时间的函数
double mytime()
{
    double ts = 0.0;
    struct timeval mt;
    gettimeofday(&mt,(struct timezone*)0);
    ts = (double)(mt.tv_sec+mt.tv_usec*1.0e-6);
    return (ts);
}


/*产生一个以0开头的 0-point_number的不重复序列*/
int *random2(int point_number)
{
	int temp = 0,signal = 1,k=0;
	int i,j;
	int *save = NULL;
	save =  (int *)malloc(point_number*sizeof(int));

	//数组初始化
	for(j = 0;j < point_number;j++)
		save[j] = -1;

	/*产生一个0-point_number的不重复序列*/
	do
	{
		signal = 1;
		temp = rand() % point_number;	//temp取0-point_number的值
		for(i=0;i<point_number;i++)
		{
			if(save[i] == temp)
			{
				signal = 0;
				break;
			}
		}
		if(signal != 0)
		{
			save[k] = temp;
			k++;
		}
	}while(signal == 0 || k != point_number);

	/*使得save[0]=0（save[0]为开始点）*/
	for(i=1;i<point_number;i++)
	{
		if(save[i]==QISHIDIAN)
		{
			temp = save[i];
			save[i] = save[0];
			save[0] = temp;
		}
	}
	return save;
}

/*产生随机数,范围在0-1之间*/
double random0_1(void)
{
	return (double)rand() / (double)RAND_MAX;
}


/*计算x y 2点之间的距离*/
double distance(int x,int y)
{
	double xy_distance = 0;
//	printf("x=%d y=%d\t",x,y);
	xy_distance = sqrt( (point_coordinate[x][0] - point_coordinate[y][0]) * (point_coordinate[x][0] - point_coordinate[y][0]) +
   						(point_coordinate[x][1] - point_coordinate[y][1]) * (point_coordinate[x][1] - point_coordinate[y][1]) );
	return xy_distance;
}


/* 产生新序列
father代表父序列，result代表子序列,point_number代表点个数 将交换结果存储在result中*/
void Neighbour(int *father,int *result,int point_number)
{
	int n = 0,m = 0,temp = 0;
	int k;
	for(k=0;k<point_number;k++)
		result[k] = father[k];

	do
	{
		n = rand() % (point_number - 1); //-1-point_number-1
		m = rand() % (point_number - 1);
	}while(n == m);		//获得随机的不相等的n，m
	n++;
	m++;//0-point_number
	//交换
	temp = result[n];
	result[n] = result[m];
	result[m] = temp;
}


/*计算某一个序列中城市之间的距离总和*/
double distance_sum(int *x,int point_number)
{
	double result = 0;
	int i;
	for(i=0;i<point_number-1;i++)
	{
		result += distance(x[i],x[i+1]);
	}
//	result += distance(x[point_number-1],x[0]);	//是否回到起始点
	return result;
}



//t代表当前温度,i,j分别代表不同的两个序列,返回对应的转移发生概率
double P(int *i,int *j,double t,int point_number)
{
	double fi = 0,fj = 0;
	double result = 0;
	fi = distance_sum(i,point_number);
	fj = distance_sum(j,point_number);
	if(fj < fi)
		result = 1.0;
	else
		result = exp( (fi - fj)/t );
	return result;
}


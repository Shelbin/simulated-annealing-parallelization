#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <sys/time.h>
#define ZUIDADIANSHU 200	//最大的点数
#define QISHIDIAN 4			//起始点
#define CHUSHIWENDU 1400000	//初始温度
#define WENDUXISHU 0.99	    //降温系数
#define T_MIN 0.000001


//全局变量
char point_list[ZUIDADIANSHU][20];				//存储点名称
double point_coordinate[ZUIDADIANSHU][2];		//存储点x y坐标

int *simulation(int *x,int point_number);
double P(int *i,int *j,double t,int point_number);
int *simulation(int *x,int point_number);
double distance_sum(int *x,int point_number);
void Neighbour(int *father,int *result,int point_number);
double distance(int x,int y);
double random0_1(void);
int *random2(int point_number);
void print_coordinate(int *list,int point_number);
double mytime();



int main()
{
//记录时间
  double tstart,tstop;

	char filename[81];
	double x_position,y_position;
	FILE *pfile;
	int *sequence;			//初始序列
	int point_number;		//点数量  可以考虑用宏定义定义最大数量
	int i;					//for循环专用变量
	int *endlist;
	double sumlen;
  tstart=mytime();

	strcpy(filename,"100.txt");
	if((pfile = fopen(filename,"r"))==NULL)
	{
		printf("没有找到城市文件%s!\n",filename);
		return 0;
	}

	/*读取目标文件城市信息*/
	fscanf(pfile,"%d",&point_number);
	for(i=0;i<point_number;i++)
	{
		fscanf(pfile,"%s\t%lf\t%lf",&point_list[i],&x_position,&y_position);
		point_coordinate[i][0] = x_position;
		point_coordinate[i][1] = y_position;
	}
	fclose(pfile);

	srand( (unsigned)time( NULL ) );   //设置种子 初始化随机发生器
	sequence = random2(point_number);

	endlist = simulation(sequence,point_number);
	printf("last: ");
	for(i=0;i<point_number;i++)
	{
		printf("%d ",endlist[i]);
	}
	sumlen = distance_sum(endlist,point_number);
	printf("\n sumlen:%lf\n",sumlen);
	free(endlist);
	free(sequence);

    tstop=mytime();
    printf("spend time=%lf\n",tstop-tstart);
	return 0;
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


/*温度的下降函数
初始温度为T=CHUSHIWENDU 温度下降系数为rate=WENDUXISHU
x为序列 point_number为点数*/
int *simulation(int *x,int point_number)
{

	int	*temp = NULL;
	double random = 0;
	int *i = malloc(point_number*sizeof(int));
	int *j = malloc(point_number*sizeof(int));
	int m;
	double t = CHUSHIWENDU, rate = WENDUXISHU;		//初始温度，降温系数

	int L = 200*point_number;			//每个温度的迭代次数,也就是每一个温度上的平衡条件
	int time = 0;						//记录某一温度下的迭代次数


	for(m=0;m<point_number;m++)
	{
		i[m] = x[m];
		j[m] = 0;
	}
double nextlen=0.0,len=0.0,local_len=0.0,d_len=0.0;
	do
	{
		do
		{
     for(m=0;m<20000;m++){
			Neighbour(i,j,point_number);
			random = P(i,j,t,point_number);
			if ( (random == 1.0) || (random > random0_1()) )
			{
				temp = i;
				i = j;
				j = temp;//保证i序列始终是当前起始序列
			 }
			}
    local_len = distance_sum(i,point_number);
    nextlen = len;
    len = local_len;
    d_len = len-nextlen;
   }while(d_len>0.1);    //		f2 = distance_sum(i,point_number);
		t = rate * t;
	}while(t>T_MIN);//结束条件2：t>t_min
	free(j);
	return i;
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
void print_coordinate(int *list,int point_number)
{
	int i;
	for(i=0;i<point_number;i++)
	{
		printf("G00 X%.3lf Y%.3lf\n",point_coordinate[list[i]][0],point_coordinate[list[i]][1]);
	}
}

double mytime()
{
    double ts = 0.0;
    struct timeval mt;
    gettimeofday(&mt,(struct timezone*)0);
    ts = (double)(mt.tv_sec+mt.tv_usec*1.0e-6);
    return (ts);
}

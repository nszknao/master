/**
 * @brief
 * 焼きなまし法
 * 円周上に100個の点を配置し，その全ての点を最短距離で通るように最適化する
 * 巡回セールスマン問題の円バージョン 
 */

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200             /* how many points do we try before stepping */
#define ITERS_FIXED_T 2000      /* how many iterations for each T? */
#define STEP_SIZE 1.0           /* max step size in random walk */
#define K 1.0                   /* Boltzmann constant */
#define T_INITIAL 5000.0        /* initial temperature */
#define MU_T 1.001              /* damping factor for temperature */
#define T_MIN 1.0e-4

/* マクロ */
#define ARRAY_LENGTH(array) (sizeof(array) / sizeof(array[0]))  // 配列の要素数を求める


gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

struct Stsp_circle {
  int number;
  double x, y;        /* coordinates */
};

void set_circle(void);
void prepare_distance_matrix(void);
void exhaustive_search(void);
void print_distance_matrix(void);
int *shuffle(int *ary, int size);
double city_distance(Stsp_circle c1, Stsp_circle c2);
double Etsp(void *xp);
double Mtsp(void *xp, void *yp);
void Stsp(const gsl_rng * r, void *xp, double step_size);
void Ptsp(void *xp);

#define N_DOTS 100    /* 円周上の点の数 */
#define RADIUS 1      /* 半径 */

struct Stsp_circle circle[N_DOTS];
double distance_matrix[N_DOTS][N_DOTS];

/**
 * @fn
 * 2点間の距離を計算
 */
double dot_distance(Stsp_circle c1, Stsp_circle c2)
{
  double distance;

  distance  = sqrt(pow((c1.x-c2.x),2) + pow((c1.y-c2.y),2));

  return distance;
}

/**
 * @fn
 * 現在の総移動距離を求める
 */
double Etsp(void *xp)
{
  int *route = (int *) xp;
  double E = 0;
  unsigned int i;

  for (i = 0; i < N_DOTS; ++i) {
    E += distance_matrix[route[i]][route[(i + 1) % N_DOTS]];
  }

  return E;
}

/**
 * @fn
 * 前後のステップを比べる
 */
double Mtsp(void *xp, void *yp)
{
  int *route1 = (int *) xp, *route2 = (int *) yp;
  double distance = 0;
  unsigned int i;

  for (i = 0; i < N_DOTS; ++i) {
    distance += ((route1[i] == route2[i]) ? 0 : 1);
  }

  return distance;
}

/**
 * @fn
 * ステップを1つすすめる
 * 配列の要素をランダムに1組だけ入れ替える
 */
void Stsp(const gsl_rng * r, void *xp, double step_size)
{
  int x1, x2, dummy;
  int *route = (int *) xp;

  step_size = 0 ; /* prevent warnings about unused parameter */

  /* 2点をランダムに選んで入れ替える */
  x1 = (gsl_rng_get (r) % (N_DOTS-1)) + 1;
  do {
    x2 = (gsl_rng_get (r) % (N_DOTS-1)) + 1;
  } while (x2 == x1);

  dummy = route[x1];
  route[x1] = route[x2];
  route[x2] = dummy;
}

/**
 * @fn
 * ステップの状況を表示する
 */
void Ptsp(void *xp)
{
  unsigned int i;
  int *route = (int *) xp;
  printf("  [");
  for (i = 0; i < N_DOTS; ++i) {
    printf(" %d ", route[i]);
  }
  printf("]  ");
}

int main(void)
{
  int *x_initial, ret[N_DOTS];
  unsigned int i;

  gsl_ieee_env_setup ();

  set_circle();
  prepare_distance_matrix();

  /* 円周上の点の初期配列を作成 */
  for (i = 0; i < N_DOTS; ++i) {
    ret[i] = i;
  }
  x_initial = shuffle(ret, N_DOTS);

  /* 各点間の距離を保持するマトリクスを表示 */
//  printf("# distance matrix is:\n");
//  print_distance_matrix();

  /* 点の初期配列を表示 */
  printf("# initial coordinates of dots\n");
  /* プロットに使うので消さない */
  for (i = 0; i < N_DOTS+1; ++i) {
    printf("###initial_dot_coord: %d %lf %lf\n",
           circle[x_initial[i % N_DOTS]].number,
           circle[x_initial[i % N_DOTS]].x,
           circle[x_initial[i % N_DOTS]].y);
  }

  const gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup()) ;
  gsl_siman_solve(r, x_initial, Etsp, Stsp, Mtsp, Ptsp, NULL, NULL, NULL,
                  N_DOTS*sizeof(int), params);

  /* 点の最終配列を表示 */
  printf("# final coordinates of dots\n");
  /* プロットに使うので消さない */
  for (i = 0; i < N_DOTS+1; ++i) {
    printf("###final_dot_coord: %d %lf %lf\n",
           circle[x_initial[i % N_DOTS]].number,
           circle[x_initial[i % N_DOTS]].x,
           circle[x_initial[i % N_DOTS]].y);
  }

  printf("# ");

  /* バッファの内容を強制的に出力（フラッシュ）する */
  fflush(stdout);
#if 0
  system("date");
#endif /* 0 */
  fflush(stdout);

  return 0;
}

/**
 * @fn
 * 円周上の各点に，番号と座標を設定する
 */
void set_circle()
{
  unsigned int i;
  double theta;

  for (int i = 0; i < N_DOTS; ++i) {
    double num_dots = (double)N_DOTS;

    theta = i/num_dots*2.*M_PI;
    circle[i].number  = i;
    circle[i].x       = RADIUS*cos(theta);
    circle[i].y       = RADIUS*sin(theta);
  }
}

/**
 * @fn
 * 各点間の距離を保持するマトリクスを作成
 */
void prepare_distance_matrix()
{
  unsigned int i, j;
  double dist;

  for (i = 0; i < N_DOTS; ++i) {
   for (j = 0; j < N_DOTS; ++j) {
      if (i == j) {
        dist = 0;
      } else {
        dist = dot_distance(circle[i], circle[j]);
      }
      distance_matrix[i][j] = dist;
    }
  }
}

/**
 * @fn
 * 各点間の距離を保持するマトリクスをコンソールに表示
 */
void print_distance_matrix()
{
  unsigned int i, j;

  for (i = 0; i < N_DOTS; ++i) {
    printf("# ");
    for (j = 0; j < N_DOTS; ++j) {
      printf("%15.8f   ", distance_matrix[i][j]);
    }
    printf("\n");
  }
}

/**
 * @fn
 * int型の配列をシャッフルする
 */
int *shuffle(int *ary, int size)
{
  unsigned int i;

    for(i=0;i<size;i++)
    {
        const gsl_rng *r = gsl_rng_alloc (gsl_rng_default);
        gsl_rng_set(r, rand());
        int j = gsl_rng_get(r)%size;
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }
    return ary;
}
#include "../include/simulation.h"
#include "../include/common.h"

/********** 計算条件 **********/
const std::size_t Simulation::SAMPLE_LENGTH = 131072; // 131072,65536
const std::size_t Simulation::NUM_OF_SAMPLES = 100; // 入力の標本数
const double Simulation::dt = 0.01; // 時間刻み幅
const std::size_t Simulation::PDF_NUM_OF_PARTITIONS = 100; // 応答分布の変位・速度分割数

Simulation::Simulation(double lambda, double beta2, double alpha):
_lambda(lambda),
_beta2(beta2),
_alpha(alpha)
{}

/**
 * @fn 入力を生成する
 * @param vector< vector<double> > &force 入力を格納（force[0]:合成した入力，force[1]:ガウス性ホワイトノイズ，force[2]:不規則パルス励振）
 **/
void Simulation::_createExcitation(std::vector< std::vector<double> > &force)
{
    std::size_t i;

    // ランダム変数
    gsl_rng* r  = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(NULL)+clock());
    gsl_rng* rp = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rp, time(NULL)+clock()+1);

    // 入力強度
    double wSt = sqrt(_alpha);
    double pSt = sqrt(1 - _alpha);

    // ホワイトノイズの分散
    double sigma = sqrt(2.*M_PI*Common::S0 / dt);

    Common::resize2DemensionalVector(force, 3, SAMPLE_LENGTH);
    #pragma omp parallel for
    for (i = 0; i < SAMPLE_LENGTH; ++i) {
        force[0][i] = wSt*gsl_ran_gaussian(r, sigma) + pSt*gsl_ran_bernoulli(r, dt*_lambda)*gsl_ran_gaussian(rp, sqrt(_beta2)) / dt;
        force[1][i] = wSt*gsl_ran_gaussian(r, sigma);
        force[2][i] = pSt*gsl_ran_bernoulli(r, dt*_lambda)*gsl_ran_gaussian(rp, sqrt(_beta2)) / dt;
    }
}

/**
 * @fn 運動方程式を4次のルンゲクッタで解く
 * @param vector<double> &x 時間を格納
 * @param vector<double> &y 応答変位を格納
 * @param vector< vector<double> > &forces 入力
 **/
void Simulation::culcRungeKutta(std::vector<double> &t, std::vector< std::vector<double> > &v_y1, std::vector< std::vector<double> > &v_y2, std::vector< std::vector<double> > &forces)
{
    std::cout << "Start solving equation of motion using 4th-Runge-Kutta.\n" << std::endl;

    std::size_t i, ii;

    t.resize(SAMPLE_LENGTH);
    Common::resize2DemensionalVector(v_y1, NUM_OF_SAMPLES, SAMPLE_LENGTH);
    Common::resize2DemensionalVector(v_y2, NUM_OF_SAMPLES, SAMPLE_LENGTH);
    // 初期値
    double y1 = 0., y2 = 0.;
    // ルンゲクッタの諸変数
    double DY1[4], DY2[4];
    for (i = 0; i < NUM_OF_SAMPLES; ++i)
    {
        // forces[0]が合成した入力
        std::vector< std::vector<double> > _forces;
        this->_createExcitation(_forces);
        std::vector<double> force;
        force   = _forces[0];

        y1 = 0., y2 = 0.;

        for (ii = 0; ii < SAMPLE_LENGTH; ++ii)
        {
            DY1[0] = dt*_f1(force[ii], y1, y2);
            DY2[0] = dt*_f2(force[ii], y1, y2);

            DY1[1] = dt*_f1(force[ii], y1 + DY1[0] / 2.0, y2 + DY2[0] / 2.0);
            DY2[1] = dt*_f2(force[ii], y1 + DY1[0] / 2.0, y2 + DY2[0] / 2.0);

            DY1[2] = dt*_f1(force[ii], y1 + DY1[1] / 2.0, y2 + DY2[1] / 2.0);
            DY2[2] = dt*_f2(force[ii], y1 + DY1[1] / 2.0, y2 + DY2[1] / 2.0);

            DY1[3] = dt*_f1(force[ii], y1 + DY1[2], y2 + DY2[2]);
            DY2[3] = dt*_f2(force[ii], y1 + DY1[2], y2 + DY2[2]);

            y1 = y1 + (DY1[0] + 2.*DY1[1] + 2.*DY1[2] + DY1[3]) / 6.0;
            y2 = y2 + (DY2[0] + 2.*DY2[1] + 2.*DY2[2] + DY2[3]) / 6.0;

            // 変位と速度の最小値・最大値を決定
            if (y1>_y1max) _y1max = y1;
            if (y1<_y1min) _y1min = y1;
            if (y2>_y2max) _y2max = y2;
            if (y2<_y2min) _y2min = y2;

            v_y1[i][ii] = y1;
            v_y2[i][ii] = y2;

            if (i == 0){
                t[ii]   = ii*dt;
                forces  = _forces;
            }
        }
    }
}

/**
 * @fn 変位の確率密度関数を生成
 * @param vector< vector<double> > &y1 応答変位が格納された２次元配列
 * @param vector<double> &x X軸の情報を保存
 * @param vector<double> &y Y軸の情報を保存
 **/
void Simulation::createDispPdf(const std::vector< std::vector<double> > &y1, std::vector<double> &x, std::vector<double> &y)
{
    std::cout << "Creating a file of the displacement pdf(.dat).\n" << std::endl;

    std::size_t i, ii, iii;

    std::size_t n_dx;
    double dx = (_y1max - _y1min)/PDF_NUM_OF_PARTITIONS;
    std::vector< std::vector<double> > pdf_buffer(NUM_OF_SAMPLES, std::vector<double>(NUM_OF_SAMPLES, PDF_NUM_OF_PARTITIONS+1));

    for (i = 0; i < NUM_OF_SAMPLES; ++i) {
        for (ii = 0; (_y1min + ii*dx) <= _y1max; ++ii) {
            // 幅dxに含まれる回数
            // TODO:バッファに追加した値を除外することで速度向上
            n_dx = 0;
            for (iii = 0; iii < SAMPLE_LENGTH; ++iii) {
                if ((y1[i][iii] < (_y1min + (ii + 1)*dx)) && (y1[i][iii] >= (_y1min + ii*dx))) n_dx++;
            }
            pdf_buffer[i][ii] = (double)n_dx / (SAMPLE_LENGTH*dx);
        }
    }

    double pdf_y1, integral_y1 = 0.;
    x.resize(PDF_NUM_OF_PARTITIONS+1);
    y.resize(PDF_NUM_OF_PARTITIONS+1);
    for (i = 0; (_y1min + i*dx) <= _y1max; ++i) {
        pdf_y1 = 0.;
        for (ii = 0; ii < NUM_OF_SAMPLES; ++ii) {
            pdf_y1 += pdf_buffer[ii][i];
        }
        pdf_y1 = (double)pdf_y1 / NUM_OF_SAMPLES;

        x[i]    = _y1min + i*dx;
        y[i]    = pdf_y1;

        integral_y1 += pdf_y1*dx;
    }

    std::cout << "integral of pdf_y1 = " << integral_y1 << ".\n" << std::endl;
}

/**
 * @fn 速度の確率密度関数を生成
 * @param vector< vector<double> > &y2 応答速度が格納された２次元配列
 * @param vector<double> &x X軸の情報を保存
 * @param vector<double> &y Y軸の情報を保存
 **/
void Simulation::createVelPdf(const std::vector< std::vector<double> > &y2, std::vector<double> &x, std::vector<double> &y)
{
    std::cout << "Creating a file of the velocity pdf(.dat).\n" << std::endl;

    std::size_t i, ii, iii;

    std::size_t n_dx;
    double dx = (_y2max - _y2min)/(PDF_NUM_OF_PARTITIONS+100);
    std::vector< std::vector<double> > pdf_buffer(NUM_OF_SAMPLES, std::vector<double>(NUM_OF_SAMPLES, PDF_NUM_OF_PARTITIONS+101));

    for (i = 0; i < NUM_OF_SAMPLES; ++i) {
        for (ii = 0; (_y2min + ii*dx) <= _y2max; ++ii) {
            // 幅dxに含まれる回数
            n_dx = 0;
            for (iii = 0; iii < SAMPLE_LENGTH; ++iii) {
                if ((y2[i][iii] < (_y2min + (ii + 1)*dx)) && (y2[i][iii] >= (_y2min + ii*dx))) n_dx++;
            }
            pdf_buffer[i][ii] = (double)n_dx / (SAMPLE_LENGTH*dx);
        }
    }

    double pdf_y2, integral_y2 = 0.;
    x.resize(PDF_NUM_OF_PARTITIONS+101);
    y.resize(PDF_NUM_OF_PARTITIONS+101);
    for (i = 0; (_y2min + i*dx) <= _y2max; ++i) {
        pdf_y2 = 0.;
        for (ii = 0; ii < NUM_OF_SAMPLES; ++ii) {
            pdf_y2 += pdf_buffer[ii][i];
        }
        pdf_y2 = (double)pdf_y2 / NUM_OF_SAMPLES;

        x[i] = _y2min + i*dx;
        y[i] = pdf_y2;

        integral_y2 += pdf_y2*dx;
    }

    std::cout << "integral of pdf_y2 = " << integral_y2 << ".\n" << std::endl;
}

/**
 * @fn 結合確率密度関数を生成
 * @param vector< vector<double> > &y1 応答変位が格納された２次元配列
 * @param vector< vector<double> > &y2 応答速度が格納された２次元配列
 * @param vector<double> &x X軸の情報を保存
 * @param vector<double> &y Y軸の情報を保存
 * @param vector<double> &z Z軸の情報を保存
 **/
void Simulation::createJointPdf(const std::vector< std::vector<double> > &y1, const std::vector< std::vector<double> > &y2, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
{
    std::cout << "Creating a file of the joint pdf(.dat).\n" << std::endl;

    std::size_t i, ii, iii, iiii;

    /* 逐次処理のスレッドに関する情報を表示 */
    std::cout << "現在使用中のスレッド数は「" << omp_get_num_threads() << "」です。" << std::endl;
    std::cout << "使用可能なスレッド数は最大「" << omp_get_max_threads() << "」です。" << std::endl;
    
    std::size_t n_dx;
    double dx_y1 = (_y1max - _y1min)/PDF_NUM_OF_PARTITIONS;
    double dx_y2 = (_y2max - _y2min)/PDF_NUM_OF_PARTITIONS;
    std::vector< std::vector< std::vector<double> > > pdf_buffer;
    Common::resize3DemensionalVector(pdf_buffer, PDF_NUM_OF_PARTITIONS+1, PDF_NUM_OF_PARTITIONS+1, NUM_OF_SAMPLES);

    for (i = 0; i<NUM_OF_SAMPLES; ++i) {
        std::cout << "sample: " << i << std::endl;
        for (ii = 0; (_y1min + ii*dx_y1) <= _y1max; ++ii) {
            for (iii = 0; (_y2min + iii*dx_y2) <= _y2max; ++iii) {
                // 幅dxに含まれる回数
                n_dx = 0;
                for (iiii = 0; iiii<SAMPLE_LENGTH; ++iiii) {
                    if ((y1[i][iiii] < (_y1min + (ii + 1)*dx_y1)) && (y1[i][iiii] >= (_y1min + ii*dx_y1)) && (y2[i][iiii] < (_y2min + (iii + 1)*dx_y2)) && (y2[i][iiii] >= (_y2min + iii*dx_y2))) n_dx++;
                }
    
                pdf_buffer[ii][iii][i] = (double)n_dx / SAMPLE_LENGTH / (dx_y1*dx_y2);
            }
        }
   }

    x.resize(PDF_NUM_OF_PARTITIONS+1);
    y.resize(PDF_NUM_OF_PARTITIONS+1);
    z.resize((PDF_NUM_OF_PARTITIONS+1)*(PDF_NUM_OF_PARTITIONS+1));
    double pdf_joint, integral_joint = 0.;
    for (i = 0; (_y1min + i*dx_y1) <= _y1max; ++i) {
        for (ii = 0; (_y2min + ii*dx_y2) <= _y2max; ++ii) {
            pdf_joint = 0.;
            for (iii = 0; iii<NUM_OF_SAMPLES; ++iii) {
                pdf_joint += pdf_buffer[i][ii][iii];
            }
            if (i == 0) {
                y[ii] = _y2min + ii*dx_y2;
            }
            pdf_joint = (double)pdf_joint / NUM_OF_SAMPLES;
            z[i*(PDF_NUM_OF_PARTITIONS+1) + ii] = pdf_joint;            
            integral_joint += pdf_joint*(dx_y1*dx_y2);
        }
        x[i] = _y1min + i*dx_y1;
    }

    std::cout << "integral of pdf_joint = " << integral_joint << ".\n" << std::endl;
}

/**
 * @fn ガウス性ホワイトノイズを受ける系の厳密解
 * @param vector<double> &x X軸の情報を保存
 * @param vector<double> &y Y軸の情報を保存
 */
void Simulation::exactSolutionOfGaussianWhiteNoise(std::vector<double> &x, std::vector<double> &y)
{
    std::cout << "Culculate exact solution excited Gaussian white noise." << std::endl;

    std::size_t i;

    // 第二修正ベッセル関数
    double K_nu;
    double bessel_p, bessel_q;
    double bessel_C, bessel_arg;

    // 確率密度関数の厳密解
    double exact_gauss_pdf  = 0.;
    double integral_gauss_pdf   = 0.;
    
    bessel_p = 2.*Common::ZETA;
    bessel_q = Common::ZETA*Common::EPSILON;
    bessel_arg = pow(bessel_p, 2) / (8.*bessel_q);

    K_nu = gsl_sf_bessel_Knu(0.25, bessel_arg);
    bessel_C = 2.*sqrt(bessel_q / bessel_p)*exp(-bessel_arg) / K_nu;
    
    double y1min = -5., y1max = 5.;
    double dx_y1 = (y1max - y1min)/PDF_NUM_OF_PARTITIONS;
    x.resize(PDF_NUM_OF_PARTITIONS);
    y.resize(PDF_NUM_OF_PARTITIONS);
    for (i = 0; (y1min + i*dx_y1) <= y1max; i++) {
        exact_gauss_pdf = bessel_C*exp(-bessel_p*pow(y1min + i*dx_y1, 2) - bessel_q*pow(y1min + i*dx_y1, 4));
        x[i] = y1min + i*dx_y1;
        y[i] = exact_gauss_pdf;
        integral_gauss_pdf += exact_gauss_pdf*dx_y1;
    }

    std::cout << "integral of exact_gauss_pdf = " << integral_gauss_pdf << ".\n" << std::endl;
}

double Simulation::_f1(double force, double y1, double y2)
{
    return y2;
}

double Simulation::_f2(double force, double y1, double y2)
{
    return force - 2.*Common::ZETA*y2 - y1 - Common::EPSILON*y1*y1*y1;
}

/**
 * main
 */
int main(int argc, char *argv[])
{
    std::string filename    = "";

    char *ends;
    double lambda   = strtod(argv[1],&ends);
    double beta2    = strtod(argv[2],&ends);
    double alpha    = strtod(argv[3],&ends);

    std::cout << "--------------------\n" << std::endl;
    std::cout << "research.cpp started.\n" << std::endl;

    Simulation * sim    = new Simulation(lambda, beta2, alpha);

    /* ルンゲクッタを解く */
    std::vector<double> t;
    std::vector< std::vector<double> > y1, y2, forces;
    sim->culcRungeKutta(t, y1, y2, forces);
    filename = "sim_force.dat";
    Common::outputIntoFile(filename, t, forces[0]);
    filename = "sim_force_gaussian.dat";
    Common::outputIntoFile(filename, t, forces[1]);
    filename = "sim_force_pulse.dat";
    Common::outputIntoFile(filename, t, forces[2]);
    filename = "sim_x1.dat";
    Common::outputIntoFile(filename, t, y1[0]);
    filename = "sim_x2.dat";
    Common::outputIntoFile(filename, t, y2[0]);
    /* 変位のPDFを求める */
    std::vector<double> dispX, dispY;
    sim->createDispPdf(y1, dispX, dispY);
    filename = "sim_y1pdf.dat";
    Common::outputIntoFile(filename, dispX, dispY);
    /* 速度のPDFを求める */
    std::vector<double> velX, velY;
    sim->createVelPdf(y2, velX, velY);
    filename = "sim_y2pdf.dat";
    Common::outputIntoFile(filename, velX, velY);
    /* 結合PDFを求める */
    std::vector<double> jointX, jointY, jointZ;
    sim->createJointPdf(y1, y2, jointX, jointY, jointZ);
    filename = "sim_jointpdf.dat";
    Common::output3DIntoFile(filename, jointX, jointY, jointZ);
    /* ホワイトノイズを受ける系の応答分布の厳密解を求める */
//    std::vector<double> gwnDispX, gwnDispY;
//    sim->exactSolutionOfGaussianWhiteNoise(gwnDispX, gwnDispY);
//    filename = "sim_y1Gpdf.dat";
//    Common::outputIntoFile(filename, gwnDispX, gwnDispY);
    delete sim;

    std::cout << "analysis.cpp has done.\n" << std::endl;
    std::cout << "--------------------\n" << std::endl;

    return 0;
}

#include "class.cuh"


void MPS::RdDat() {
	//////////////////////physical.txt///////////////////////////////
	FILE* in;
	if (fopen_s(&in, IN_FILE_1, "r") != 0) {
		printf_s("MPS_physical.txtが開けません\n");
	}
	else {
		real scan[50];
		fscanf_s(in, "%f %f %f", &scan[0], &scan[1], &scan[2]);//最小範囲
		fscanf_s(in, "%f %f %f", &scan[3], &scan[4], &scan[5]);//最大範囲
		fscanf_s(in, "%f %f %f", &scan[6], &scan[7], &scan[8]);//重力加速度
		fscanf_s(in, "%f", &scan[9]);//計算終了時間
		fscanf_s(in, "%f", &scan[10]);//ファイル吐き出し感覚
		fscanf_s(in, "%f", &scan[11]);//マッハ数
		fscanf_s(in, "%f", &scan[12]);//クーラン数
		fscanf_s(in, "%f", &scan[13]);//剛体衝突の反発係数
		fscanf_s(in, "%f", &scan[14]);//限界接近距離
		fscanf_s(in, "%f", &scan[15]);//流体の動粘性係数
		fscanf_s(in, "%f", &scan[16]);//粒子サイズ
		fclose(in);

		PCL_DST = scan[16];
		MINc.x = real(scan[0] - 3.1 * PCL_DST);
		MINc.y = real(scan[1] - 3.1 * PCL_DST);
		MINc.z = real(scan[2] - 3.1 * PCL_DST);
		MAXc.x = real(scan[3] + 3.1 * PCL_DST);
		MAXc.y = real(scan[4] + 3.1 * PCL_DST);
		MAXc.z = real(scan[5] + 3.1 * PCL_DST);
		G.x = real(scan[6]);
		G.y = real(scan[7]);
		G.z = real(scan[8]);
		FIN_TIM = real(scan[9]);
		output_time = real(scan[10]);
		Ma = real(scan[11]);
		CRT_NUM = real(scan[12]);
		COL_RAT = real(scan[13]);
		DST_LMT_RAT = real(scan[14]);
		KNM_VSC = real(scan[15]);
	}
	fclose(in);

	printf_s("MINc.x = %f  MINc.y = %f  MINc.z = %f\n", MINc.x, MINc.y, MINc.z);
	printf_s("MAXc.x = %f  MAXc.y = %f  MAXc.z = %f\n", MAXc.x, MAXc.y, MAXc.z);
	printf_s("G.x = %f  G.y = %f  G.z = %f\n", G.x, G.y, G.z);
	printf_s("FIN_TIM = %f\n", FIN_TIM);
	printf_s("output_time = %f\n", output_time);
	printf_s("Ma = %f\n", Ma);
	printf_s("CRT_NUM = %f\n", CRT_NUM);
	printf_s("COL_RAT = %f\n", COL_RAT);
	printf_s("DST_LMT_RAT = %f\n", DST_LMT_RAT);
	printf_s("KNM_VSC = %f\n", KNM_VSC);
	printf_s("PCL_DST = %f\n", PCL_DST);
	//////////////////////physical.txt///////////////////////////////

	//////////////////////initial.txt///////////////////////////////
	FILE* in2;
	if (fopen_s(&in2, IN_FILE_2, "r") != 0) {
		printf_s("initial.txtが開けません\n");
	}
	else {
		fscanf_s(in2, "%d", &nP);//総粒子数取得

		//流体
		Pos.x = (real*)malloc(sizeof(real) * (nP));
		Pos.y = (real*)malloc(sizeof(real) * (nP));
		Pos.z = (real*)malloc(sizeof(real) * (nP));

		Vel.x = (real*)malloc(sizeof(real) * (nP));
		Vel.y = (real*)malloc(sizeof(real) * (nP));
		Vel.z = (real*)malloc(sizeof(real) * (nP));

		Acc.x = (real*)malloc(sizeof(real) * (nP));
		Acc.y = (real*)malloc(sizeof(real) * (nP));
		Acc.z = (real*)malloc(sizeof(real) * (nP));

		Prs = (real*)malloc(sizeof(real) * (nP));
		pav = (real*)malloc(sizeof(real) * (nP));
		//n = (real*)malloc(sizeof(real) * (nP));
		//流体


		Typ = (char*)malloc(sizeof(char) * (nP));
		Dns = (real*)malloc(sizeof(real) * (Dns_Num));

		//壁
		WLLVec.x = (real*)malloc(sizeof(real) * (nP));
		WLLVec.y = (real*)malloc(sizeof(real) * (nP));
		WLLVec.z = (real*)malloc(sizeof(real) * (nP));
		WLLSE = (char*)malloc(sizeof(char) * (nP));
		//壁

		//ミラー
		PosM.x = (real*)malloc(sizeof(real) * (nP * NumMRR));
		PosM.y = (real*)malloc(sizeof(real) * (nP * NumMRR));
		PosM.z = (real*)malloc(sizeof(real) * (nP * NumMRR));
		VelM.x = (real*)malloc(sizeof(real) * (nP * NumMRR));
		VelM.y = (real*)malloc(sizeof(real) * (nP * NumMRR));
		VelM.z = (real*)malloc(sizeof(real) * (nP * NumMRR));
		PrsM = (real*)malloc(sizeof(real) * (nP * NumMRR));
		TypM = (char*)malloc(sizeof(char) * (nP * NumMRR));


		int nPfluid = 0;
		int nPwall = 0;
		int nPobj = 0;
		int nPtmp = 0;
		for (int i = 0; i < nP; i++) {
			int a[1];
			float b[11];
			int c[1];
			float g[1];
			fscanf_s(in2, " %d %d %f %f %f %f %f %f %f", &a[0], &c[0], &b[0], &b[1], &b[2], &b[8], &b[9], &b[10], &g[0]);
			const treal3 pos = { b[0], b[1], b[2] };
			if (pos.x<MAXc.x && pos.x>MINc.x && pos.y<MAXc.y && pos.y>MINc.y && pos.z<MAXc.z && pos.z>MINc.z) {
				Typ[nPtmp] = char(c[0]);
				Pos.x[nPtmp] = real(b[0]); Pos.y[nPtmp] = real(b[1]); Pos.z[nPtmp] = real(b[2]);
				Vel.x[nPtmp] = Vel.y[nPtmp] = Vel.z[nPtmp] = 0.0f;
				Acc.x[nPtmp] = Acc.y[nPtmp] = Acc.z[nPtmp] = 0.0f;
				Prs[nPtmp] = 0.0f;
				pav[nPtmp] = 0.0f;
				//n[nPtmp] = 0.0f;

				if (Typ[nPtmp] == FLD) { nPfluid += 1; WLLVec.x[nPtmp] = 0.0f; WLLVec.y[nPtmp] = 0.0f; WLLVec.z[nPtmp] = 0.0f; }
				else if (Typ[nPtmp] == WLL) { nPwall += 1; WLLVec.x[nPtmp] = -real(b[8]); WLLVec.y[nPtmp] = -real(b[9]); WLLVec.z[nPtmp] = -real(b[10]); }
				else if (Typ[nPtmp] == OBJ) { nPobj += 1; WLLVec.x[nPtmp] = -real(b[8]); WLLVec.y[nPtmp] = -real(b[9]); WLLVec.z[nPtmp] = -real(b[10]); }
				nPtmp += 1;
			}
		}
		nP = nPtmp;
		nPWLL = nPwall;
		nPFLD = nPfluid;
		nPOBJ = nPobj;
		std::cout << "総流体粒子数 nPFLD = " << nPFLD << std::endl;
		std::cout << "総壁粒子数 nPWLL = " << nPWLL << std::endl;
		std::cout << "総動壁粒子数 nPOBJ = " << nPOBJ << std::endl;
		std::cout << "総粒子数 nP = " << nP << std::endl;
	}
	fclose(in2);
	//////////////////////initial.txt///////////////////////////////
}


__global__ void d_initialize_int_array(const int n, int* i_array, const int a) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n) {
		i_array[i] = a;
	}
}

__global__ void d_initialize_real_array(const int n, real* i_array, real a) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n) {
		i_array[i] = a;
	}
}


void MPS::Output_vtk(const char typ, const int outputflg) {
	if (outputflg == 0)
	{
		sprintf_s(outout_filename, "./vtk_cuda/output%05d.csv", iF);
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf("%sが開けません\n", outout_filename);
		}
		/*else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,Vel,Prs,pav\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%f\n", Pos.x[i], Pos.y[i], Pos.z[i], sqrt(Vel.x[i] * Vel.x[i] + Vel.y[i] * Vel.y[i] + Vel.z[i] * Vel.z[i]), Prs[i], pav[i] / OPT_FQC);
				}
			}
		}*/
		else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,Vel.x,Vel.y,Vel.z,Vel,Prs,pav\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f\n", Pos.x[i], Pos.y[i], Pos.z[i], Vel.x[i], Vel.y[i], Vel.z[i], sqrt(Vel.x[i] * Vel.x[i] + Vel.y[i] * Vel.y[i] + Vel.z[i] * Vel.z[i]), Prs[i], pav[i] / OPT_FQC);
				}
			}
		}
		fclose(fp);
	}

	if (outputflg == 1)
	{
		sprintf_s(outout_filename, "./vtk_cuda/outputwall%05d.csv", iF);
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf_s("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,Prs,pav,SE\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%d\n", Pos.x[i], Pos.y[i], Pos.z[i], Prs[i], pav[i] / OPT_FQC, WLLSE[i]);
				}
			}
		}
		fclose(fp);
	}

	if (outputflg == 2)
	{
		sprintf_s(outout_filename, "./vtk_cuda/outputMRR%05d.csv", iF);
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf_s("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "PosM.x,PosM.y,PosM.z,VelM,PrsM\n");
			for (int i = 0; i < nP; i++) {
				for (int k = 0; k < NumMRR; k++) {
					int kiNM = k + i * NumMRR;
					if (TypM[kiNM] == typ) {
						fprintf_s(fp, "%f,%f,%f,%f,%f,\n", PosM.x[kiNM], PosM.y[kiNM], PosM.z[kiNM], sqrt(VelM.x[kiNM] * VelM.x[kiNM] + VelM.y[kiNM] * VelM.y[kiNM] + VelM.z[kiNM] * VelM.z[kiNM]), PrsM[kiNM]);
					}
				}
			}
		}
		fclose(fp);
	}

	if (outputflg == 3)
	{
		sprintf_s(outout_filename, "./vtk_cuda/outputWLLvec%05d.csv", iF);
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf_s("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "PosM.x,PosM.y,PosM.z,vec.x,vec.y,vec.z\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%f\n", Pos.x[i], Pos.y[i], Pos.z[i], WLLVec.x[i], WLLVec.y[i], WLLVec.z[i]);
				}
			}
		}
		fclose(fp);
	}
}


void MPS::Output_vtk2(const char typ, const int outputflg) {//txt seitei
	if (outputflg == 0)
	{
		sprintf_s(outout_filename, "./vtk_cuda/seitei.txt");
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf_s("%sが開けません\n", outout_filename);
		}
		else {
			int fldnp = 0;
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == FLD) { fldnp += 1; }//流体粒子カウント
			}
			fprintf_s(fp, "%d\n", fldnp);
			int FLDnP = 0;
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == FLD) {
					fprintf_s(fp, "%d %d %f %f %f %f %f %f %f\n", FLDnP, Typ[i], Pos.x[i], Pos.y[i], Pos.z[i], 0.0, 0.0, 0.0, 0.0);
					FLDnP += 1;
				}
			}
		}
		fclose(fp);
	}
}


void MPS::WrtDat(void) {//流体出力
	Output_vtk(FLD, 0);//水粒子出力
	Output_vtk(MRR, 2);//生成ミラー確認
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	((d_initialize_real_array << <blocks_nP, threads >> > (nP, d_pav, 0.0f)));//平均圧力初期化
	CHECK(cudaDeviceSynchronize());
	printf_s("WrtDat_done!\n\n");
}

void MPS::WrtDatWLL(void) {//壁出力
	Output_vtk(WLL, 1);//壁出力
	Output_vtk(WLL, 3);//法線ベクトル
	printf("WrtDatWLL_done!\n\n");
}


void MPS::WrtDat2(void) {//静定状態出力
	Output_vtk2(FLD, 0);
	printf_s("WrtDat2_done!\n\n");
}


void MPS::AlcBkt() {//バケット作成

	r = PCL_DST * 3.1f;
	r2 = r * r;
	rp = PCL_DST * 2.1f;//圧力用
	rp2 = rp * rp;
	DB = PCL_DST * 3.1f;
	DB2 = DB * DB;
	DBinv = 1.0f / DB;

	nBx = (int)((MAXc.x - MINc.x) * DBinv) + 3;
	nBy = (int)((MAXc.y - MINc.y) * DBinv) + 3;
	nBz = (int)((MAXc.z - MINc.z) * DBinv) + 3;

	nBxy = nBx * nBy;
	nBxyz = nBx * nBy * nBz;
	printf_s("nBx:%d  nBy:%d  nBz:%d  nBxy:%d  nBxyz:%d\n", nBx, nBy, nBz, nBxy, nBxyz);

	(cudaMalloc((void**)&d_bfst, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_blst, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_nxt, sizeof(int) * (nP)));
	(cudaMalloc((void**)&d_bfstM, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_blstM, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_nxtM, sizeof(int) * (nP * NumMRR)));
}


void MPS::SetPara() {
	//初期粒子密度
	real tn0 = 0.0f;
	real tn0_grad = 0.0f;
	real tlmd = 0.0f;
	real tlmd_grad = 0.0f;
	int tNumn0_grad = 0;
	for (int ix = -10; ix < 10; ix++) {
		for (int iy = -10; iy < 10; iy++) {
			for (int iz = -10; iz < 10; iz++) {
				real x = real(PCL_DST) * real(ix);
				real y = real(PCL_DST) * real(iy);
				real z = real(PCL_DST) * real(iz);
				real dist2 = x * x + y * y + z * z;
				if (dist2 == 0.0f) { continue; }
				real dist = sqrt(dist2);
				if (dist2 < rp2) {
					tn0_grad += WEI_grad(dist, rp);
					tlmd_grad += dist2 * WEI_grad(dist, rp);
					tNumn0_grad += 1;
				}
				if (dist2 < r2) {
					tn0 += WEI(dist, r);
					tlmd += dist2 * WEI(dist, r);
				}
			}
		}
	}
	n0 = tn0;
	lmd = tlmd / tn0;
	n0_grad = tn0_grad;
	lmd_grad = tlmd_grad / tn0_grad;
	Numn0_grad = tNumn0_grad;
	printf_s("n0:%f\nn0_grad:%f\nNumn0_grad:%d\n", n0, n0_grad, Numn0_grad);

	//接近禁止距離
	rlim = PCL_DST * DST_LMT_RAT;
	rlim2 = rlim * rlim;

	COL = real(COL_RAT + 1.0f);//反発係数設定

#pragma omp parallel for
	for (int i = 0; i < Dns_Num; i++) {//初期化
		Dns[i] = 0.0f;
	}
	Dns[FLD] = Dns_FLD;
	Dns[WLL] = Dns_WLL;
	Dns[OBJ] = Dns_OBJ;
	Dns[OBJ2] = Dns_OBJ2;
	Dns[MRR] = Dns_FLD;


#pragma omp parallel for
	for (int i = 0; i < nP; i++) {//初期化
		for (int k = 0; k < NumMRR; k++) {
			int kiNM = k + i * NumMRR;
			PosM.x[kiNM] = PosM.y[kiNM] = PosM.z[kiNM] = 0.0f;
			VelM.x[kiNM] = VelM.y[kiNM] = VelM.z[kiNM] = 0.0f;
			PrsM[kiNM] = 0.0f;
			TypM[kiNM] = GST;
		}
		WLLSE[i] = Surface;
	}

	//音速とタイムステップ
	real max_height = MINc.y;
	real min_height = MAXc.y;
	for (int i = 0; i < nP; i++) {//初期化
		if (Typ[i] == FLD) {
			if (max_height < Pos.y[i]) { max_height = Pos.y[i]; }
			else if (min_height > Pos.y[i]) { min_height = Pos.y[i]; }
		}
	}
	umax = sqrt(2.0f * abs(G.y) * (max_height - min_height));
	SND = umax / Ma;
	Prs_coef = SND * SND / n0_grad;
	Vsc_coef = 2.0f * 3.0f * KNM_VSC / n0 / lmd;
	Pmax = 1000.0f * G.y* (max_height - min_height);//9.80665f 

	//dt = CRT_NUM * PCL_DST / SND;
	//dt = 0.000005f;//固定
	dt = real(PCL_DST / (umax + SND));//計算間隔(ダムブレイクの最大速度を基準にする)

	printf_s("タイムステップ:dt=%f\n最大速度:umax=%f\n音速=%f\nVsc_coef=%f\nPrs_coef=%f\nPmax=%f\n", dt, umax, SND, Vsc_coef, Prs_coef, Pmax);

	//配列確保とVRAMへの転送
	(cudaMalloc((void**)&d_Typ, sizeof(char) * nP));
	(cudaMalloc((void**)&d_Pos.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Pos.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Pos.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Vel.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Vel.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Vel.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Acc.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Acc.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Acc.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Prs, sizeof(real) * nP));
	(cudaMalloc((void**)&d_pav, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Dns, sizeof(real) * Dns_Num));

	(cudaMalloc((void**)&d_WLLVec.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_WLLVec.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_WLLVec.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_WLLSE, sizeof(char) * nP));

	(cudaMalloc((void**)&d_FromWLL, sizeof(int) * (nP * NumMRR)));//デバイスのみ(転送なし)

	(cudaMalloc((void**)&d_TypM, sizeof(char) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PosM.x, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PosM.y, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PosM.z, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_VelM.x, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_VelM.y, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_VelM.z, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PrsM, sizeof(real) * (nP * NumMRR)));

	(cudaMemcpy(d_Typ, Typ, sizeof(char) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Pos.x, Pos.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Pos.y, Pos.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Pos.z, Pos.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Vel.x, Vel.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Vel.y, Vel.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Vel.z, Vel.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Acc.x, Acc.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Acc.y, Acc.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Acc.z, Acc.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Prs, Prs, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_pav, pav, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Dns, Dns, sizeof(real) * Dns_Num, cudaMemcpyHostToDevice));

	(cudaMemcpy(d_WLLVec.x, WLLVec.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_WLLVec.y, WLLVec.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_WLLVec.z, WLLVec.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_WLLSE, WLLSE, sizeof(char) * nP, cudaMemcpyHostToDevice));

	(cudaMemcpy(d_TypM, TypM, sizeof(char) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PosM.x, PosM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PosM.y, PosM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PosM.z, PosM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_VelM.x, VelM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_VelM.y, VelM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_VelM.z, VelM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PrsM, PrsM, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));

}


__device__ void d_ChkPcl(const int i, char* d_Typ, areal3 d_Pos, areal3 d_Vel, const  areal3 d_Acc, real* d_Prs, const  treal3 MINc, const  treal3 MAXc, const real PCL_DST, const real umax)//&要る？
{
	if (d_Typ[i] == FLD) {//最大速度に制限
		if (d_Pos.x[i] < (MAXc.x - 3.1f * PCL_DST) && d_Pos.x[i] > (MINc.x + 3.1f * PCL_DST) &&
			d_Pos.y[i] < (MAXc.y - 3.1f * PCL_DST) && d_Pos.y[i] > (MINc.y + 3.1f * PCL_DST) &&
			d_Pos.z[i] < (MAXc.z - 3.1f * PCL_DST) && d_Pos.z[i] > (MINc.z + 3.1f * PCL_DST)) {


			treal3 Utmp;
			Utmp.x = d_Vel.x[i];	Utmp.y = d_Vel.y[i];	Utmp.z = d_Vel.z[i];
			real U = Utmp.x * Utmp.x + Utmp.y * Utmp.y + Utmp.z * Utmp.z;
			U = sqrt(U);
			if (U > umax) {
				Utmp.x *= umax / U;	Utmp.y *= umax / U;	Utmp.z *= umax / U;
				d_Vel.x[i] = Utmp.x;		d_Vel.y[i] = Utmp.y;		d_Vel.z[i] = Utmp.z;
			}


			/*real u_limit = 1.5f * umax;
			if (d_Vel.x[i] > u_limit) { d_Vel.x[i] = u_limit; }
			if (d_Vel.y[i] > u_limit) { d_Vel.y[i] = u_limit; }
			if (d_Vel.z[i] > u_limit) { d_Vel.z[i] = u_limit; }*/

		}
		else {//発散している粒子は削除
			d_Typ[i] = GST;
			d_Pos.x[i] = d_Pos.y[i] = d_Pos.z[i] = 0.0f;
			d_Vel.x[i] = d_Vel.y[i] = d_Vel.z[i] = 0.0f;
			d_Acc.x[i] = d_Acc.y[i] = d_Acc.z[i] = 0.0f;
			d_Prs[i] = 0.0f;
		}
	}
}


__global__ void d_MkBkt(const int nP, const  int nBx, const  int nBxy, const  real DBinv, int* d_bfst, int* d_blst, int* d_nxt, const char* d_Typ, areal3 d_Pos, const treal3 MINc)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		if (d_Typ[i] == GST) { return; }
		int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
		int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
		int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
		int ib = iz * nBxy + iy * nBx + ix;
		const int j = atomicExch(&d_blst[ib], i);
		if (j == -1) { d_bfst[ib] = i; }
		else { d_nxt[j] = i; }
	}
}


void MPS::MkBkt() {//粒子をバケットに収納
	//printf_s("MkBkt start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_bfst, -1)));
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_blst, -1)));
	((d_initialize_int_array << <blocks_nP, threads >> > (nP, d_nxt, -1)));
	CHECK(cudaDeviceSynchronize());

	d_MkBkt << <blocks_nP, threads >> > (nP, nBx, nBxy, DBinv, d_bfst, d_blst, d_nxt, d_Typ, d_Pos, MINc);
	CHECK(cudaDeviceSynchronize());

	//printf_s("MkBkt done!\n\n");
}


__global__ void d_Surface_Edge(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_WLLVec, char* d_WLLSE, const real PCL_DST,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)//壁粒子surface edge判定
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		real r2 = (1.0f * PCL_DST) * (1.0f * PCL_DST);//隣接壁粒子判定距離の2乗
		if (d_Typ[i] == WLL) {
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			treal3 vec;
			vec.x = d_WLLVec.x[i];	vec.y = d_WLLVec.y[i];	vec.z = d_WLLVec.z[i];

			int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
			int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
			int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == WLL) {
									treal3 p;//i,jの距離の成分
									p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < r2) {
										if ((d_WLLVec.x[j] != vec.x) || (d_WLLVec.y[j] != vec.y) || (d_WLLVec.z[j] != vec.z)) {//近くの粒子と法線ベクトル違うならedge
											d_WLLSE[i] = Edge;
										}
									}
								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
		}

		if (d_Typ[i] == OBJ) {
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			treal3 vec;
			vec.x = d_WLLVec.x[i];	vec.y = d_WLLVec.y[i];	vec.z = d_WLLVec.z[i];

			int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
			int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
			int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == OBJ) {
									treal3 p;//i,jの距離の成分
									p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < r2) {
										if ((d_WLLVec.x[j] != vec.x) || (d_WLLVec.y[j] != vec.y) || (d_WLLVec.z[j] != vec.z)) {//近くの粒子と法線ベクトル違うならedge
											d_WLLSE[i] = Edge;
										}
									}
								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
		}

		if (d_Typ[i] == OBJ2) {
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			treal3 vec;
			vec.x = d_WLLVec.x[i];	vec.y = d_WLLVec.y[i];	vec.z = d_WLLVec.z[i];

			int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
			int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
			int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == OBJ2) {
									treal3 p;//i,jの距離の成分
									p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < r2) {
										if ((d_WLLVec.x[j] != vec.x) || (d_WLLVec.y[j] != vec.y) || (d_WLLVec.z[j] != vec.z)) {//近くの粒子と法線ベクトル違うならedge
											d_WLLSE[i] = Edge;
										}
									}
								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
		}

	}
}


void MPS::Surface_Edge() {//surface-edge判定
	//printf_s("Surface_Edge  start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_Surface_Edge << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_WLLVec, d_WLLSE, PCL_DST, MINc, DBinv, nBx, nBxy, d_bfst, d_blst, d_nxt);
	CHECK(cudaDeviceSynchronize());

	CHECK(cudaMemcpy(WLLSE, d_WLLSE, sizeof(char) * nP, cudaMemcpyDeviceToHost));

	//printf_s("Surface_Edge done!\n\n");
}




__global__ void d_ResetMRR(const int nP_NumMRR, char* d_TypM, areal3 d_PosM, areal3 d_VelM, real* d_PrsM, int* d_FromWLL)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP_NumMRR) {
		d_TypM[i] = GST;
		d_PosM.x[i] = d_PosM.y[i] = d_PosM.z[i] = 0.0f;
		d_VelM.x[i] = d_VelM.y[i] = d_VelM.z[i] = 0.0f;
		d_PrsM[i] = 0.0f;
		d_FromWLL[i] = -1;
	}
}


void MPS::ResetMRR() {//ミラー粒子削除
	//printf_s("ResetMRR start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP * NumMRR);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP_NumMRR(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_ResetMRR << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), d_TypM, d_PosM, d_VelM, d_PrsM, d_FromWLL);
	CHECK(cudaDeviceSynchronize());

	//printf_s("ResetMRR done!\n\n");
}


__global__ void d_GenMRR_nonslip(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, real* d_Prs, areal3 d_WLLVec, char* d_TypM, areal3 d_PosM, areal3 d_VelM, real* d_PrsM, int* d_FromWLL, char* d_WLLSE, const real PCL_DST, const real r2,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		char Typ = d_Typ[i];
		if (Typ == FLD) {//流体粒子から面対称ミラー生成
			real r2_2 = 1.0f * r2;
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			int WLLexist = 0;//近傍に壁粒子が存在したら１
			int iNM = i * NumMRR;
			int Edgeexist = 0;
			int Edge_unique[NumMRR];
			for (int k = 0; k < NumMRR; k++) { Edge_unique[k] = -1; }

			//最近傍壁粒子探索
			int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
			int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
			int iz = (int)((pos.z - MINc.z) * DBinv) + 1;
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								char typ = d_Typ[j];
								if ((typ == WLL) || (typ == OBJ) || (typ == OBJ2)) {
									treal3 p;
									p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];//壁粒子→流体粒子のベクトル
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
									if (dist2 <= r2_2) {//影響半径内に壁有り
										WLLexist = 1;//壁有りフラグ
										treal3 wallvec;
										wallvec.x = d_WLLVec.x[j];	wallvec.y = d_WLLVec.y[j];	wallvec.z = d_WLLVec.z[j];

										int unique = 0;//面のミラー生成
										for (int k = 0; k < NumMRR; k++) {//unique判定(iからk番目に生成されたミラーと法線ベクトルが重複していないか判定)
											int jNum = d_FromWLL[k + iNM];
											if (jNum != -1) {
												if ((d_WLLVec.x[jNum] == wallvec.x) && (d_WLLVec.y[jNum] == wallvec.y) && (d_WLLVec.z[jNum] == wallvec.z)) {
													unique = 1;//同じ法線ベクトルが登録済み													
												}
											}
										}
										if (unique == 0) {//jがuniqueであれば登録
											for (int k = 0; k < NumMRR; k++) {
												int kiNM = k + iNM;
												if (d_FromWLL[kiNM] == -1) {
													d_FromWLL[kiNM] = j;
													break;
												}
											}
										}

										if (d_WLLSE[j] == Edge) {//角のミラー生成
											Edgeexist = 1;
											int E_unique = 0;
											for (int k = 0; k < NumMRR; k++) {
												int q = Edge_unique[k];
												if (q != -1) {
													if ((d_WLLVec.x[q] == wallvec.x) && (d_WLLVec.y[q] == wallvec.y) && (d_WLLVec.z[q] == wallvec.z)) {
														E_unique = 1;//同じ法線ベクトルが登録済み														
													}
												}
											}
											if (E_unique == 0) {
												for (int k = 0; k < NumMRR; k++) {
													if (Edge_unique[k] == -1) {
														Edge_unique[k] = j;
														break;
													}
												}
											}
										}

									}
								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}

			//ミラー粒子生成
			if (WLLexist == 1) {//近くに壁がある

				for (int k = 0; k < NumMRR; k++) {
					int kiNM = k + iNM;
					int FromNum = d_FromWLL[kiNM];//壁粒子番号jをFromNumとしてレジスタに登録
					if (FromNum == -1) { continue; }//壁粒子j(FromNum)からミラー生成	
					treal3 posw;
					posw.x = d_Pos.x[FromNum];	posw.y = d_Pos.y[FromNum];	posw.z = d_Pos.z[FromNum];
					treal3 WLLvec;
					WLLvec.x = d_WLLVec.x[FromNum];	WLLvec.y = d_WLLVec.y[FromNum];	WLLvec.z = d_WLLVec.z[FromNum];
					treal3 PWvec;//最近傍壁粒子と流体粒子iの相対座標ベクトル
					PWvec.x = pos.x - posw.x;		PWvec.y = pos.y - posw.y;		PWvec.z = pos.z - posw.z;

					real PW_WLL = PWvec.x * WLLvec.x + PWvec.y * WLLvec.y + PWvec.z * WLLvec.z;
					if (PW_WLL < 0) { continue; }//内外判定　内側にはミラー生成しない(流体が壁の外にある)

					real distance = PW_WLL / sqrt(WLLvec.x * WLLvec.x + WLLvec.y * WLLvec.y + WLLvec.z * WLLvec.z);//距離のスカラー値 法線ベクトルなので分母は１になってるはず
					d_TypM[kiNM] = MRR;
					d_PosM.x[kiNM] = posw.x + PWvec.x - 2.0f * distance * WLLvec.x;//法線ベクトルを用いて対称位置に生成
					d_PosM.y[kiNM] = posw.y + PWvec.y - 2.0f * distance * WLLvec.y;
					d_PosM.z[kiNM] = posw.z + PWvec.z - 2.0f * distance * WLLvec.z;
					/*d_VelM.x[kiNM] = -d_Vel.x[i];//non-slip
					d_VelM.y[kiNM] = -d_Vel.y[i];
					d_VelM.z[kiNM] = -d_Vel.z[i];*/
					d_PrsM[kiNM] = d_Prs[i];
					d_VelM.x[kiNM] =0.0f;//non-slip
					d_VelM.y[kiNM] =0.0f;
					d_VelM.z[kiNM] =0.0f;
					//d_PrsM[kiNM] = 0.0f;
				}//ここまで面ミラー生成

#if 1 //Edge生成フラグ
				if (Edgeexist == 1) {
					int NumEdge = 0;
					for (int k = 0; k < NumMRR; k++) {
						if (Edge_unique[k] != -1) { NumEdge++; }
					}

					if (NumEdge >= 2) {//Edge粒子が影響半径内に2種類以上あったら合成ベクトルから角のミラー生成
						treal3 posw;
						treal3 synvec;
						treal3 PWvec;
						real distance;
						for (int k = 0; k < NumMRR; k++) {//合成ベクトル計算+最後のedgeを生成元に設定	
							int FromNum = Edge_unique[k];
							if (FromNum == -1) { continue; }
							posw.x = d_Pos.x[FromNum];	posw.y = d_Pos.y[FromNum];	posw.z = d_Pos.z[FromNum];
							synvec.x += d_WLLVec.x[FromNum];		synvec.y += d_WLLVec.y[FromNum];		synvec.z += d_WLLVec.z[FromNum];
							PWvec.x = pos.x - posw.x;		PWvec.y = pos.y - posw.y;		PWvec.z = pos.z - posw.z;
						}

						real abs_synvec = sqrt(synvec.x * synvec.x + synvec.y * synvec.y + synvec.z * synvec.z);
						synvec.x /= abs_synvec;	synvec.y /= abs_synvec;	synvec.z /= abs_synvec;//合成ベクトルの大きさを１にする

						real PW_WLL = PWvec.x * synvec.x + PWvec.y * synvec.y + PWvec.z * synvec.z;
						if (PW_WLL < 0) { return; }//内外判定　内側にはミラー生成しない(流体が壁の外にある)

						distance = PW_WLL / sqrt(synvec.x * synvec.x + synvec.y * synvec.y + synvec.z * synvec.z);//距離のスカラー値 法線ベクトルなので分母は１になってるはず

						int MRR_space = -1;
						for (int k = 0; k < NumMRR; k++) {
							int kiNM = k + iNM;
							if (d_FromWLL[kiNM] == -1) { MRR_space = kiNM; break; }
						}

						d_TypM[MRR_space] = MRR;
						d_PosM.x[MRR_space] = posw.x + PWvec.x - 2.0f * distance * synvec.x;//法線ベクトルを用いて対称位置に生成
						d_PosM.y[MRR_space] = posw.y + PWvec.y - 2.0f * distance * synvec.y;
						d_PosM.z[MRR_space] = posw.z + PWvec.z - 2.0f * distance * synvec.z;
						/*d_VelM.x[MRR_space] = -d_Vel.x[i];
						d_VelM.y[MRR_space] = -d_Vel.y[i];
						d_VelM.z[MRR_space] = -d_Vel.z[i];
						d_PrsM[MRR_space] = d_Prs[i];*/
						d_VelM.x[MRR_space] = 0.0f;
						d_VelM.y[MRR_space] = 0.0f;
						d_VelM.z[MRR_space] = 0.0f;
						d_PrsM[MRR_space] = d_Prs[i];
					}
				}
#endif
			}

		}
	}
}


void MPS::GenMRR_nonslip() {//ミラー全成分生成
	//printf_s("GenMRR_nonslip start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP * NumMRR);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nPMRR(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_GenMRR_nonslip << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Prs, d_WLLVec, d_TypM, d_PosM, d_VelM, d_PrsM, d_FromWLL, d_WLLSE, PCL_DST, r2, MINc, DBinv, nBx, nBxy, d_bfst, d_blst, d_nxt);
	CHECK(cudaDeviceSynchronize());

	//printf_s("GenMRR_nonslip done!\n\n");
}




__global__ void d_MkBkt_MRR(const int nP_NumMRR, const  int nBx, const  int nBxy, const  real DBinv, int* d_bfstM, int* d_blstM, int* d_nxtM, const char* d_TypM, const  areal3 d_PosM, const treal3 MINc)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP_NumMRR) {
		if (d_TypM[i] == GST) { return; }
		int ix = (int)((d_PosM.x[i] - MINc.x) * DBinv) + 1;
		int iy = (int)((d_PosM.y[i] - MINc.y) * DBinv) + 1;
		int iz = (int)((d_PosM.z[i] - MINc.z) * DBinv) + 1;
		int ib = iz * nBxy + iy * nBx + ix;
		const int j = atomicExch(&d_blstM[ib], i);
		if (j == -1) { d_bfstM[ib] = i; }
		else { d_nxtM[j] = i; }
	}
}

//壁の内側にできたミラーをはじく
__global__ void d_MRRinout(const int nP_NumMRR, const  int nBx, const  int nBxy, const  real DBinv, int* d_bfst, int* d_blst, int* d_nxt, int* d_bfstM, int* d_blstM, int* d_nxtM, char* d_TypM, const  areal3 d_PosM, const char* d_Typ, const areal3 d_Pos, const areal3 d_WLLVec, const real r, const treal3 MINc)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP_NumMRR) { return; }

	treal3 posM;
	posM.x = d_PosM.x[i];	posM.y = d_PosM.y[i];	posM.z = d_PosM.z[i];
	int WLLexist = 0;
	real walldist = 10.0f * r;
	treal3 posw = { 0.0f };
	treal3 wallvec = { 0.0f };

	int ix = (int)((posM.x - MINc.x) * DBinv) + 1;
	int iy = (int)((posM.y - MINc.y) * DBinv) + 1;
	int iz = (int)((posM.z - MINc.z) * DBinv) + 1;

	//最近傍壁粒子探索
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始

					char typ = d_Typ[j];
					if ((typ == WLL) || (typ == OBJ) || (typ == OBJ2)) {
						treal3 posj;
						posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];
						treal3 p;
						p.x = posM.x - posj.x;	p.y = posM.y - posj.y;	p.z = posM.z - posj.z;
						real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
						real dist = sqrt(dist2);
						if (dist <= r) {//影響半径内に壁有り
							if (dist < walldist) {
								WLLexist = 1;//壁有りフラグ jにする
								walldist = dist;
								wallvec.x = d_WLLVec.x[j];	wallvec.y = d_WLLVec.y[j];	wallvec.z = d_WLLVec.z[j];
								posw.x = posj.x;	posw.y = posj.y;	posw.z = posj.z;
							}
						}
					}

					j = d_nxt[j];
					if (j == -1) break;
				}
			}
		}
	}

	real inout = (posM.x - posw.x) * wallvec.x + (posM.y - posw.y) * wallvec.y + (posM.z - posw.z) * wallvec.z;//壁粒子→ミラー粒子のベクトル * 壁法線ベクトル
	if ((WLLexist == 1) && (inout >= 0)) { d_TypM[i] = GST; }//内外判定に引っかかったらGSTに

}

void MPS::MkBkt_MRR() {//粒子をバケットに収納
	//printf_s("MkBkt_MRR start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP * NumMRR);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP_NumMRR(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_bfstM, -1)));
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_blstM, -1)));
	((d_initialize_int_array << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), d_nxtM, -1)));
	d_MkBkt_MRR << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), nBx, nBxy, DBinv, d_bfstM, d_blstM, d_nxtM, d_TypM, d_PosM, MINc);

	d_MRRinout << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), nBx, nBxy, DBinv, d_bfst, d_blst, d_nxt, d_bfstM, d_blstM, d_nxtM, d_TypM, d_PosM, d_Typ, d_Pos, d_WLLVec, r, MINc);

	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_bfstM, -1)));
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_blstM, -1)));
	((d_initialize_int_array << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), d_nxtM, -1)));
	d_MkBkt_MRR << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), nBx, nBxy, DBinv, d_bfstM, d_blstM, d_nxtM, d_TypM, d_PosM, MINc);
	CHECK(cudaDeviceSynchronize());

	//printf_s("MkBkt_MRR done!\n\n");
}


__global__ void d_VscTrm(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, areal3 d_WLLVec, areal3 d_PosM, areal3 d_VelM, const real r, const real PCL_DST, const real n0, const real KNM_VSC, const real Vsc_coef, treal3 G,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		if (d_Typ[i] == FLD) {
			real vsc_coef = Vsc_coef;
			treal3 pos;		treal3 vel;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			vel.x = d_Vel.x[i];		vel.y = d_Vel.y[i];		vel.z = d_Vel.z[i];
			treal3 Acc;
			Acc.x = Acc.y = Acc.z = 0.0f;//加速度の一時計算
			int WLLexist = 0;
			real walldist = 10.0f * r;
			treal3 posw = { 0.0f };
			treal3 wallvec = { 0.0f };
			real invn0 = 1.0f / n0;

			int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
			int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
			int iz = (int)((pos.z - MINc.z) * DBinv) + 1;



			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2;// i,j距離の2乗
									dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
									real dist;//粒子の距離(絶対値)
									dist = sqrt(dist2);
									if (dist <= r) {
										/*real w = WEI(dist, r);//i,jの重み
										Acc.x += (d_Vel.x[j] - vel.x) * w;
										Acc.y += (d_Vel.y[j] - vel.y) * w;
										Acc.z += (d_Vel.z[j] - vel.z) * w;*/
										Acc.x += 2.0f * (d_Vel.x[j] - vel.x) * PCL_DST / dist / dist / dist;
										Acc.y += 2.0f * (d_Vel.y[j] - vel.y) * PCL_DST / dist / dist / dist;
										Acc.z += 2.0f * (d_Vel.z[j] - vel.z) * PCL_DST / dist / dist / dist;
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
			//ミラーループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfstM[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							//内外判定
							treal3 posMj;
							posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];


							treal3 p;//i,jの距離の成分
							p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
							real dist = sqrt(dist2);//粒子の距離(絶対値)
							if (dist <= r) {
								/*real w = WEI(dist, r);//i,jの重み
								Acc.x += (d_VelM.x[j] - vel.x) * w;
								Acc.y += (d_VelM.y[j] - vel.y) * w;
								Acc.z += (d_VelM.z[j] - vel.z) * w;*/
								Acc.x += 2.0f * (d_VelM.x[j] - vel.x) * PCL_DST / dist / dist / dist;
								Acc.y += 2.0f * (d_VelM.y[j] - vel.y) * PCL_DST / dist / dist / dist;
								Acc.z += 2.0f * (d_VelM.z[j] - vel.z) * PCL_DST / dist / dist / dist;
							}

							j = d_nxtM[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}

			/*d_Acc.x[i] = vsc_coef * Acc.x + G.x;
			d_Acc.y[i] = vsc_coef * Acc.y + G.y;
			d_Acc.z[i] = vsc_coef * Acc.z + G.z;*/

			real n0KNM = invn0 * KNM_VSC;
			d_Acc.x[i] = n0KNM * Acc.x + G.x;
			d_Acc.y[i] = n0KNM * Acc.y + G.y;
			d_Acc.z[i] = n0KNM * Acc.z + G.z;
		}
	}
}


void MPS::VscTrm() {//粘性項・外力項の計算
	//printf_s("VscTrm start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_VscTrm << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_WLLVec, d_PosM, d_VelM, r, PCL_DST, n0, KNM_VSC, Vsc_coef, G, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	CHECK(cudaDeviceSynchronize());

	//printf_s("VscTrm done!\n\n");
}


__global__ void d_UpPcl1(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, real* d_Prs, const  treal3 MINc, const  treal3 MAXc, const  real dt, const real  PCL_DST, const real umax)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		if (d_Typ[i] == FLD) {
			d_Vel.x[i] += d_Acc.x[i] * dt;	d_Vel.y[i] += d_Acc.y[i] * dt;	d_Vel.z[i] += d_Acc.z[i] * dt;
			d_Pos.x[i] += d_Vel.x[i] * dt;	d_Pos.y[i] += d_Vel.y[i] * dt;	d_Pos.z[i] += d_Vel.z[i] * dt;
			d_Acc.x[i] = d_Acc.y[i] = d_Acc.z[i] = 0.0f;
		}
		d_ChkPcl(i, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, MINc, MAXc, PCL_DST, umax);
	}
}


void MPS::UpPcl1() {//仮の粒子移動
	//printf_s("UpPcl1 start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_UpPcl1 << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, MINc, MAXc, dt, PCL_DST, umax);
	CHECK(cudaDeviceSynchronize());

	//printf_s("UpPcl1 done!\n\n");
}


__global__ void d_ChkCol(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, areal3 d_WLLVec, char* d_TypM, areal3 d_PosM, areal3 d_VelM, real* d_Dns, const real PCL_DST, const real r, const real rlim2, const real COL,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		char Typ = d_Typ[i];
		if (Typ == FLD) {
			treal3 pos;		treal3 vec;	treal3 vec2;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			vec.x = d_Vel.x[i];		vec.y = d_Vel.y[i];		vec.z = d_Vel.z[i];
			vec2.x = d_Vel.x[i];		vec2.y = d_Vel.y[i];		vec2.z = d_Vel.z[i];
			real mi = d_Dns[Typ];
			int WLLexist = 0;
			real walldist = 10.0f * r;
			treal3 posw = { 0.0f };
			treal3 wallvec = { 0.0f };
			real rlim_wall = 0.45f * PCL_DST;
			real rlim_wall2 = rlim_wall * rlim_wall;

			int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
			int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
			int iz = (int)((pos.z - MINc.z) * DBinv) + 1;



			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								char typ = d_Typ[j];
								if (typ == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < rlim2) {
										real fDT = (vec.x - d_Vel.x[j]) * p.x + (vec.y - d_Vel.y[j]) * p.y + (vec.z - d_Vel.z[j]) * p.z;
										if (fDT > 0.0f) {
											real mj = d_Dns[typ];
											fDT *= COL * mj / (mi + mj) / dist2;
											vec2.x -= p.x * fDT;	vec2.y -= p.y * fDT;	vec2.z -= p.z * fDT;
										}
									}

								}
								else if ((typ == WLL) || (typ == OBJ) || (typ == OBJ2)) {

									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];
									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < rlim_wall2) {
										real fDT = (vec.x - d_Vel.x[j]) * p.x + (vec.y - d_Vel.y[j]) * p.y + (vec.z - d_Vel.z[j]) * p.z;
										if (fDT > 0.0f) {
											real mj = d_Dns[typ];
											fDT *= COL * mj / (mi + mj) / dist2;
											vec2.x -= p.x * fDT;	vec2.y -= p.y * fDT;	vec2.z -= p.z * fDT;
										}
									}
								}

							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
			//ミラーループ
			real mj = d_Dns[MRR];
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfstM[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							//内外判定
							treal3 posMj;
							posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];

							treal3 p;//i,jの距離の成分
							p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
							if (dist2 < rlim2) {
								real fDT = (vec.x - d_VelM.x[j]) * p.x + (vec.y - d_VelM.y[j]) * p.y + (vec.z - d_VelM.z[j]) * p.z;
								if (fDT > 0.0f) {
									fDT *= COL * mj / (mi + mj) / dist2;
									vec2.x -= p.x * fDT;	vec2.y -= p.y * fDT;	vec2.z -= p.z * fDT;
								}
							}

							j = d_nxtM[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
			d_Acc.x[i] = vec2.x;	d_Acc.y[i] = vec2.y;	d_Acc.z[i] = vec2.z;
		}
	}
}


void __global__ d_UpChkCol(const int nP, char* d_Typ, areal3 d_Vel, areal3 d_Acc) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		if (d_Typ[i] != FLD) { return; }
		d_Vel.x[i] = d_Acc.x[i];		d_Vel.y[i] = d_Acc.y[i];		d_Vel.z[i] = d_Acc.z[i];
	}
}

void MPS::ChkCol() {//仮の粒子移動
	//printf_s("ChkCol start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_ChkCol << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_WLLVec, d_TypM, d_PosM, d_VelM, d_Dns, PCL_DST, r, rlim2, COL, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	d_UpChkCol << <blocks_nP, threads >> > (nP, d_Typ, d_Vel, d_Acc);
	CHECK(cudaDeviceSynchronize());

	//printf_s("ChkCol done!\n\n");
}


__global__ void d_MkPrs(const int nP, char* d_Typ, areal3 d_Pos, real* d_Dns, real* d_Prs, areal3 d_WLLVec, char* d_TypM, areal3 d_PosM, const real rp, const real n0_grad, const int Numn0_grad, const real Pmax, const real Prs_coef,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		char Typ = d_Typ[i];
		if (Typ == FLD) {
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			real ni = 0.0f;
			int tNumn0_grad = 0;
			int WLLexist = 0;
			real walldist = 10.0f * rp;
			treal3 posw = { 0.0f };
			treal3 wallvec = { 0.0f };

			int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
			int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
			int iz = (int)((pos.z - MINc.z) * DBinv) + 1;



			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									real dist = sqrt(dist2);;
									if (dist < rp) {
										ni += WEI_grad(dist, rp);
										tNumn0_grad += 1;
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}

			//ミラーループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfstM[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							//内外判定
							treal3 posMj;
							posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];

							treal3 p;//i,jの距離の成分
							p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
							real dist = sqrt(dist2);
							if (dist < rp) {
								ni += WEI_grad(dist, rp);
								tNumn0_grad += 1;
							}

							j = d_nxtM[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}

			real mi = d_Dns[Typ];
			real pressure = (ni > n0_grad) * (ni - n0_grad) * Prs_coef * mi;
			//real pressure = (tNumn0_grad > 0.8f*Numn0_grad) * (ni - n0_grad) * Prs_coef * mi;//個数判定

			if (pressure <= Pmax && pressure >= 0.0f)
			{
				//	pressure = Pmax;//正常値なら何もしない
			}
			else
			{
				if (pressure > Pmax)
				{
					pressure = Pmax;//最大値抑制
				}
				else
				{
					pressure = 0.0f;//負圧を0に
				}
			}

			d_Prs[i] = pressure;
		}
	}
}


void MPS::MkPrs() {//仮の粒子移動
	//printf_s("MkPrs start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_MkPrs << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Dns, d_Prs, d_WLLVec, d_TypM, d_PosM, rp, n0_grad, Numn0_grad, Pmax, Prs_coef, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	CHECK(cudaDeviceSynchronize());

	//printf_s("MkPrs done!\n\n");
}

/// <summary>
/// 入部仲座
///  成分全部足し合わせてからインバースじゃない
/// 各ループごとにインバースまで計算し得た幸せ？
/// </summary>
#if 0
__global__ void d_PrsGrdTrm(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Acc, real* d_Prs, areal3 d_WLLVec, real* d_Dns, areal3 d_PosM, real* d_PrsM, const real rp, const real rp2, const real n0_grad,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		char Typ = d_Typ[i];
		if (Typ == FLD) {
			real invn0 = 1.0f / n0_grad;
			real invro = -1.0f / d_Dns[Typ];//-1/ρ係数
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			real Prs_min = d_Prs[i];
			/*real invLeft[9] = {0.0f};//圧力勾配項左ブロックのインバースする前の3*3テンソル
			real Left[9] = { 0.0f };//圧力勾配左ブロックの3*3テンソル
			real Right[3] = { 0.0f };//圧力勾配右ブロックの圧力かけてるベクトル*/

			real invLeft0 = 0.0f;
			real invLeft1 = 0.0f;
			real invLeft2 = 0.0f;
			real invLeft3 = 0.0f;
			real invLeft4 = 0.0f;
			real invLeft5 = 0.0f;
			real invLeft6 = 0.0f;
			real invLeft7 = 0.0f;
			real invLeft8 = 0.0f;

			real Left0 = 0.0f;
			real Left1 = 0.0f;
			real Left2 = 0.0f;
			real Left3 = 0.0f;
			real Left4 = 0.0f;
			real Left5 = 0.0f;
			real Left6 = 0.0f;
			real Left7 = 0.0f;
			real Left8 = 0.0f;

			real Right0 = 0.0f;
			real Right1 = 0.0f;
			real Right2 = 0.0f;

			int WLLexist = 0;
			real walldist2 = (10.0f * rp) * (10.0f * rp);
			treal3 posw = { 0.0f };
			treal3 wallvec = { 0.0f };

			int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
			int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
			int iz = (int)((pos.z - MINc.z) * DBinv) + 1;


			//近傍最小圧力抽出
			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < rp2) {
										real prs = d_Prs[j];
										if (prs < Prs_min) {
											Prs_min = prs;
										}
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}


			//圧力項計算  入部・仲座モデル
			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2;// i,j距離の2乗
									dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
									if (dist2 < rp2) {
										real dist = sqrt(dist2);
										real w = WEI_grad(dist, rp) / dist2;
										real Prsj_min = d_Prs[j] - Prs_min;
										//real Prsj_min = d_Prs[j] + d_Prs[i];
										invLeft0 += w * p.x * p.x;
										invLeft1 += w * p.x * p.y;
										invLeft2 += w * p.x * p.z;
										invLeft3 += w * p.y * p.x;
										invLeft4 += w * p.y * p.y;
										invLeft5 += w * p.y * p.z;
										invLeft6 += w * p.z * p.x;
										invLeft7 += w * p.z * p.y;
										invLeft8 += w * p.z * p.z;
										Right0 += w * Prsj_min * p.x;
										Right1 += w * Prsj_min * p.y;
										Right2 += w * Prsj_min * p.z;
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
			//ミラーループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfstM[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								//内外判定
								treal3 posMj;
								posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];


								treal3 p;//i,jの距離の成分
								p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
								real dist2;// i,j距離の2乗
								dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
								if (dist2 < rp2) {
									real dist = sqrt(dist2);
									real w = WEI_grad(dist, rp) / dist2;
									real PrsMj_min = d_PrsM[j] - Prs_min;
									//real PrsMj_min = d_PrsM[j] + d_Prs[i];
									invLeft0 += w * p.x * p.x;
									invLeft1 += w * p.x * p.y;
									invLeft2 += w * p.x * p.z;
									invLeft3 += w * p.y * p.x;
									invLeft4 += w * p.y * p.y;
									invLeft5 += w * p.y * p.z;
									invLeft6 += w * p.z * p.x;
									invLeft7 += w * p.z * p.y;
									invLeft8 += w * p.z * p.z;
									Right0 += w * PrsMj_min * p.x;
									Right1 += w * PrsMj_min * p.y;
									Right2 += w * PrsMj_min * p.z;
								}

							}
							j = d_nxtM[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}


			invLeft0 = invn0 * invLeft0;//n0で割る
			invLeft1 = invn0 * invLeft1;
			invLeft2 = invn0 * invLeft2;
			invLeft3 = invn0 * invLeft3;
			invLeft4 = invn0 * invLeft4;
			invLeft5 = invn0 * invLeft5;
			invLeft6 = invn0 * invLeft6;
			invLeft7 = invn0 * invLeft7;
			invLeft8 = invn0 * invLeft8;

			if ((invLeft0 == 0.0f) && (invLeft1 == 0.0f) && (invLeft2 == 0.0f) && (invLeft3 == 0.0f) && (invLeft4 == 0.0f) && (invLeft5 == 0.0f) && (invLeft6 == 0.0f) && (invLeft7 == 0.0f) && (invLeft8 == 0.0f)) { return; }//ディターミナント計算で0割しないように

			real DET_Left = 1.0f / ((invLeft0 * invLeft4 * invLeft8 + invLeft1 * invLeft5 * invLeft6 + invLeft2 * invLeft3 * invLeft7) - (invLeft2 * invLeft4 * invLeft6 + invLeft0 * invLeft5 * invLeft7 + invLeft1 * invLeft3 * invLeft8));//ここからインバース計算
			Left0 = DET_Left * (invLeft4 * invLeft8 - invLeft5 * invLeft7);
			Left1 = DET_Left * -(invLeft1 * invLeft8 - invLeft2 * invLeft7);
			Left2 = DET_Left * (invLeft1 * invLeft5 - invLeft2 * invLeft4);
			Left3 = DET_Left * -(invLeft3 * invLeft8 - invLeft5 * invLeft6);
			Left4 = DET_Left * (invLeft0 * invLeft8 - invLeft2 * invLeft6);
			Left5 = DET_Left * -(invLeft0 * invLeft5 - invLeft2 * invLeft3);
			Left6 = DET_Left * (invLeft3 * invLeft7 - invLeft4 * invLeft6);
			Left7 = DET_Left * -(invLeft0 * invLeft7 - invLeft1 * invLeft6);
			Left8 = DET_Left * (invLeft0 * invLeft4 - invLeft1 * invLeft3);

			Right0 = invn0 * Right0;//n0で割る
			Right1 = invn0 * Right1;
			Right2 = invn0 * Right2;



			d_Acc.x[i] = invro * (Left0 * Right0 + Left1 * Right1 + Left2 * Right2);
			d_Acc.y[i] = invro * (Left3 * Right0 + Left4 * Right1 + Left5 * Right2);
			d_Acc.z[i] = invro * (Left6 * Right0 + Left7 * Right1 + Left8 * Right2);


		}
	}
}

#else
__global__ void d_PrsGrdTrm(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Acc, real* d_Prs, areal3 d_WLLVec, real* d_Dns, areal3 d_PosM, real* d_PrsM, const real rp, const real rp2, const real n0_grad,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		char Typ = d_Typ[i];
		if (Typ == FLD) {
			real invn0 = 1.0f / n0_grad;
			real invro = -1.0f / d_Dns[Typ];//-1/ρ係数
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			real Prs_min = d_Prs[i];

			treal3 Acc = { 0.0f };
			real A3 = 3.0f / n0_grad;//Dimension / n0

			int WLLexist = 0;
			real walldist2 = (10.0f * rp) * (10.0f * rp);
			treal3 posw = { 0.0f };
			treal3 wallvec = { 0.0f };

			int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
			int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
			int iz = (int)((pos.z - MINc.z) * DBinv) + 1;



			//近傍最小圧力抽出
			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < rp2) {
										real prs = d_Prs[j];
										if (prs < Prs_min) {
											Prs_min = prs;
										}
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}

			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {

									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2;// i,j距離の2乗
									dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
									if (dist2 < rp2) {
										real dist = sqrt(dist2);
										real w = WEI_grad(dist, rp) / dist2;
										real Prsj_min = d_Prs[j] - Prs_min;
										Acc.x += Prsj_min * w * p.x;
										Acc.y += Prsj_min * w * p.y;
										Acc.z += Prsj_min * w * p.z;
										/*real Prs = d_Prs[j] +d_Prs[i];
										Acc.x += Prs * w * p.x;
										Acc.y += Prs * w * p.y;
										Acc.z += Prs * w * p.z;*/
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
			//ミラーループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfstM[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								//内外判定
								treal3 posMj;
								posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];

								treal3 p;//i,jの距離の成分
								p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
								real dist2;// i,j距離の2乗
								dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
								if (dist2 < rp2) {
									real dist = sqrt(dist2);
									real w = WEI_grad(dist, rp) / dist2;
									real PrsMj_min = d_PrsM[j] - Prs_min;
									Acc.x += PrsMj_min * w * p.x;
									Acc.y += PrsMj_min * w * p.y;
									Acc.z += PrsMj_min * w * p.z;
									/*real Prs = d_PrsM[j] + d_Prs[i];
									Acc.x += Prs * w * p.x;
									Acc.y += Prs * w * p.y;
									Acc.z += Prs * w * p.z;*/
								}

							}
							j = d_nxtM[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}

			d_Acc.x[i] = invro * Acc.x * A3;
			d_Acc.y[i] = invro * Acc.y * A3;
			d_Acc.z[i] = invro * Acc.z * A3;
		}
	}
}
#endif


void MPS::PrsGrdTrm() {//仮の粒子移動

	//printf_s("PrsGrdTrm start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_PrsGrdTrm << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Acc, d_Prs, d_WLLVec, d_Dns, d_PosM, d_PrsM, rp, rp2, n0_grad, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	CHECK(cudaDeviceSynchronize());

	//printf_s("PrsGrdTrm done!\n\n");
}


__global__ void d_UpPcl2(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, real* d_Prs, real* d_pav, const  treal3 MINc, const  treal3 MAXc, const  real dt, const real  PCL_DST, const real umax)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		if (d_Typ[i] == FLD) {
			treal3 acc;
			acc.x = d_Acc.x[i];	acc.y = d_Acc.y[i];	acc.z = d_Acc.z[i];
			d_Vel.x[i] += acc.x * dt;	d_Vel.y[i] += acc.y * dt;	d_Vel.z[i] += acc.z * dt;
			d_Pos.x[i] += acc.x * dt * dt;	d_Pos.y[i] += acc.y * dt * dt;	d_Pos.z[i] += acc.z * dt * dt;
			d_Acc.x[i] = d_Acc.y[i] = d_Acc.z[i] = 0.0f;
		}
		d_pav[i] += d_Prs[i];//時間平均圧力加算(壁も更新)
		d_ChkPcl(i, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, MINc, MAXc, PCL_DST, umax);
	}
}

void MPS::UpPcl2() {//圧力修正分の粒子移動
	//printf_s("UpPcl2 start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_UpPcl2 << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, d_pav, MINc, MAXc, dt, PCL_DST, umax);
	CHECK(cudaDeviceSynchronize());

	//printf_s("UpPcl2 done!\n\n");
}


void MPS::DevicetoHost() {
	//printf_s("DevicetoHost start!\n");
	CHECK(cudaMemcpy(Typ, d_Typ, sizeof(char) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Pos.x, d_Pos.x, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Pos.y, d_Pos.y, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Pos.z, d_Pos.z, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Vel.x, d_Vel.x, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Vel.y, d_Vel.y, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Vel.z, d_Vel.z, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Prs, d_Prs, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(pav, d_pav, sizeof(real) * nP, cudaMemcpyDeviceToHost));

	CHECK(cudaMemcpy(TypM, d_TypM, sizeof(char) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PosM.x, d_PosM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PosM.y, d_PosM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PosM.z, d_PosM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(VelM.x, d_VelM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(VelM.y, d_VelM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(VelM.z, d_VelM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PrsM, d_PrsM, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));

	//printf_s("DevicetoHost done!\n\n");
}


void MPS::ClcMPS() {
	iF = 0;
	TIM = 0.0f;
	outtime = 0.0f;
	OPT_FQC = 1.0f;

	MkBkt();
	Surface_Edge();//壁動かすときは毎ステップ実行できるように移動させること！

	WrtDatWLL();
	WrtDat();
	printf_s("\n");
	iF++;
	OPT_FQC = 0.0f;
	while (1) {

		if (outtime >= output_time) {
			outtime -= output_time;
			//ColForce2();//壁の圧力求める？
			DevicetoHost();

			WrtDat();
			//printf_s("outtime!!!\n\n");
			printf_s("Time = %f\n\n", TIM);
			iF++;
			if (TIM >= FIN_TIM) {
				WrtDat2();
				break;
			}
			OPT_FQC = 0.0f;
		}

		MkBkt();

		ResetMRR();
		GenMRR_nonslip();
		MkBkt_MRR();

		/*evicetoHost();
			WrtDat();////////////////////////////////////////////奇数
			iF++;*/
		VscTrm();

		UpPcl1();

		MkBkt();

		ChkCol();

		MkBkt();//いらん？

		MkPrs();

		ResetMRR();
		GenMRR_nonslip();
		MkBkt_MRR();

		/*DevicetoHost();
		WrtDat();////////////////////////////////////////////偶数
		iF++;*/

		PrsGrdTrm();

		UpPcl2();

		//printf_s("TIM = %f\nouttime = %f\n\n", TIM, outtime);

		outtime += dt;
		TIM += dt;
		OPT_FQC += 1.0f;

	}

}
#pragma once

#include "const.cuh"

class MPS {
public:
	///////////関数////////////////
	void RdDat();
	void AlcBkt();
	void SetPara();
	void ClcMPS();
	void MkBkt();
	void VscTrm();
	void UpPcl1();
	void ChkCol();
	void MkPrs();
	void PrsGrdTrm();
	void UpPcl2();
	void WrtDat();
	void WrtDatWLL();
	void WrtDat2();
	void Output_vtk(char, int);
	void Output_vtk2(char, int);
	void DevicetoHost();
	void Surface_Edge();
	void ResetMRR();
	void GenMRR_nonslip();
	void MkBkt_MRR();
	///////////関数////////////////


	///////////physical.txt////////////////
	treal3 MINc;//計算範囲min
	treal3 MAXc;//計算範囲max
	treal3 G;//重力
	real FIN_TIM;//計算終了タイム
	real output_time;//ファイル出力時間間隔
	real Ma;//マッハ数
	real CRT_NUM;//クーラン数
	real COL_RAT;//剛体衝突の反発係数
	real COL;
	real DST_LMT_RAT;//限界接近距離
	real KNM_VSC;//流体の動粘性係数
	real PCL_DST;//粒子サイズ
	///////////physical.tst////////////////


	///////////粒子データ////////////////
	areal3 Pos;//座標
	areal3 d_Pos;
	areal3 Vel;//速度
	areal3 d_Vel;
	areal3 Acc;//加速度
	areal3 d_Acc;
	real* Prs;//圧力
	real* d_Prs;
	real* pav;//圧力平均
	real* d_pav;
	real* n;//粒子数密度
	real* d_n;
	real* Dns;//粒子密度
	real* d_Dns;

	char* Typ;//粒子タイプ　固体０　壁１
	char* d_Typ;

	char* WLLSE;
	char* d_WLLSE;//壁タイプ　surface　edge
	areal3 WLLVec;//壁の法線ベクトル(内向き)
	areal3 d_WLLVec;


	//ミラー粒子
	areal3 PosM;//座標
	areal3 d_PosM;
	areal3 VelM;//速度
	areal3 d_VelM;
	real* PrsM;//圧力
	real* d_PrsM;

	char* TypM;//ミラータイプ　ミラー3　ゴースト-１
	char* d_TypM;

	int* d_FromWLL;//どの壁粒子(番号)から生成されたミラーかを保存


	/// ////////////////////////////////////////////////////////////////////////////


	real r;//影響半径
	real r2;//影響半径の2乗
	real rp;
	real rp2;
	real rlim;//接近禁止半径
	real rlim2;//接近禁止半径の2乗

	real n0;//初期粒子密度
	real lmd;//ラプラシアン係数
	real n0_grad;
	real lmd_grad;

	int Numn0_grad;

	real umax;//限界速度
	real SND;//音速
	real Prs_coef;//仮の圧力を求めるときの係数
	real Vsc_coef;//粘性項の係数
	real Pmax;

	int nP;//初期粒子数
	int nPWLL;
	int nPFLD;
	int nPOBJ;
	int nPOBJ2;
	int nPMRR;

	real w;//回転壁角速度
	real* wallangle;//壁回転角
	real* d_wallangle;
	real* rot_rad;//壁回転半径
	real* d_rot_rad;
	///////////粒子データ////////////////


	///////////出力////////////////
	real dt;//タイムステップ
	int iF = 0;//ファイル番号
	real TIM = 0.0f;//現在の計算時間
	real outtime = 0.0f;
	char outout_filename[256];
	FILE* fp;
	real OPT_FQC = 1.0f;
	///////////出力////////////////


	///////////粒子探索用////////////////
	real DB, DB2, DBinv;//粒子探索の分割法のブロックの幅
	int nBx, nBy, nBz, nBxy, nBxyz;//リンクリスト用
	//int* bfst, * blst, * nxt;
	int* d_bfst, * d_blst, * d_nxt;//リンクリスト用
	int* d_bfstM, * d_blstM, * d_nxtM;//ミラーリンクリスト用
	///////////粒子探索用////////////////

};
#pragma once

#include "const.cuh"

class MPS {
public:
	///////////�֐�////////////////
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
	///////////�֐�////////////////


	///////////physical.txt////////////////
	treal3 MINc;//�v�Z�͈�min
	treal3 MAXc;//�v�Z�͈�max
	treal3 G;//�d��
	real FIN_TIM;//�v�Z�I���^�C��
	real output_time;//�t�@�C���o�͎��ԊԊu
	real Ma;//�}�b�n��
	real CRT_NUM;//�N�[������
	real COL_RAT;//���̏Փ˂̔����W��
	real COL;
	real DST_LMT_RAT;//���E�ڋߋ���
	real KNM_VSC;//���̂̓��S���W��
	real PCL_DST;//���q�T�C�Y
	///////////physical.tst////////////////


	///////////���q�f�[�^////////////////
	areal3 Pos;//���W
	areal3 d_Pos;
	areal3 Vel;//���x
	areal3 d_Vel;
	areal3 Acc;//�����x
	areal3 d_Acc;
	real* Prs;//����
	real* d_Prs;
	real* pav;//���͕���
	real* d_pav;
	real* n;//���q�����x
	real* d_n;
	real* Dns;//���q���x
	real* d_Dns;

	char* Typ;//���q�^�C�v�@�ő̂O�@�ǂP
	char* d_Typ;

	char* WLLSE;
	char* d_WLLSE;//�ǃ^�C�v�@surface�@edge
	areal3 WLLVec;//�ǂ̖@���x�N�g��(������)
	areal3 d_WLLVec;


	//�~���[���q
	areal3 PosM;//���W
	areal3 d_PosM;
	areal3 VelM;//���x
	areal3 d_VelM;
	real* PrsM;//����
	real* d_PrsM;

	char* TypM;//�~���[�^�C�v�@�~���[3�@�S�[�X�g-�P
	char* d_TypM;

	int* d_FromWLL;//�ǂ̕Ǘ��q(�ԍ�)���琶�����ꂽ�~���[����ۑ�


	/// ////////////////////////////////////////////////////////////////////////////


	real r;//�e�����a
	real r2;//�e�����a��2��
	real rp;
	real rp2;
	real rlim;//�ڋߋ֎~���a
	real rlim2;//�ڋߋ֎~���a��2��

	real n0;//�������q���x
	real lmd;//���v���V�A���W��
	real n0_grad;
	real lmd_grad;

	int Numn0_grad;

	real umax;//���E���x
	real SND;//����
	real Prs_coef;//���̈��͂����߂�Ƃ��̌W��
	real Vsc_coef;//�S�����̌W��
	real Pmax;

	int nP;//�������q��
	int nPWLL;
	int nPFLD;
	int nPOBJ;
	int nPOBJ2;
	int nPMRR;

	real w;//��]�Ǌp���x
	real* wallangle;//�ǉ�]�p
	real* d_wallangle;
	real* rot_rad;//�ǉ�]���a
	real* d_rot_rad;
	///////////���q�f�[�^////////////////


	///////////�o��////////////////
	real dt;//�^�C���X�e�b�v
	int iF = 0;//�t�@�C���ԍ�
	real TIM = 0.0f;//���݂̌v�Z����
	real outtime = 0.0f;
	char outout_filename[256];
	FILE* fp;
	real OPT_FQC = 1.0f;
	///////////�o��////////////////


	///////////���q�T���p////////////////
	real DB, DB2, DBinv;//���q�T���̕����@�̃u���b�N�̕�
	int nBx, nBy, nBz, nBxy, nBxyz;//�����N���X�g�p
	//int* bfst, * blst, * nxt;
	int* d_bfst, * d_blst, * d_nxt;//�����N���X�g�p
	int* d_bfstM, * d_blstM, * d_nxtM;//�~���[�����N���X�g�p
	///////////���q�T���p////////////////

};
#include "class.cuh"

int main(int argc, char** argv) {
	std::chrono::system_clock::time_point  start, end; // �^�� auto �ŉ�
	start = std::chrono::system_clock::now(); // �v���J�n����

	MPS obj;

	printf("�g�p�\�ȍő�X���b�h���F%d\n", omp_get_max_threads());
#pragma omp parallel for
	for (int i = 0; i < omp_get_max_threads(); i++) {
		printf("start MPS_omp_%d\n", i);
	}

	printf("RdDat\n\n");
	obj.RdDat();

	printf("AlkBkt\n\n");
	obj.AlcBkt();

	printf("SetPara\n\n");
	obj.SetPara();

	obj.ClcMPS();

	printf("end MPS\n\n");

	end = std::chrono::system_clock::now();  // �v���I������
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); //�����ɗv�������Ԃ�b�ɕϊ�
	std::cout << elapsed << " second" << std::endl;

	int t = elapsed;
	int h = t / 3600;   t %= 3600;
	int m = t / 60;     t %= 60;
	int s = t;
	std::cout << h << "h " << m << "m " << s << "s " << std::endl;

	return 0;
}
#include "class.cuh"

int main(int argc, char** argv) {
	std::chrono::system_clock::time_point  start, end; // 型は auto で可
	start = std::chrono::system_clock::now(); // 計測開始時間

	MPS obj;

	printf("使用可能な最大スレッド数：%d\n", omp_get_max_threads());
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

	end = std::chrono::system_clock::now();  // 計測終了時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); //処理に要した時間を秒に変換
	std::cout << elapsed << " second" << std::endl;

	int t = elapsed;
	int h = t / 3600;   t %= 3600;
	int m = t / 60;     t %= 60;
	int s = t;
	std::cout << h << "h " << m << "m " << s << "s " << std::endl;

	return 0;
}
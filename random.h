#include <ctime>
#include <random>

std::ranlux24_base gen(1455115);
std::uniform_real_distribution<> dist1(-0.5,0.5);
std::uniform_real_distribution<> dist(0.0, 1.0);

double rand_zero_one()
{
	return dist(gen);
}
double rand_half()
{
	return dist1(gen);
}
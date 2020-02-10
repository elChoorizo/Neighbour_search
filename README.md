# Neighbour_search
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
using namespace std;


struct Particle {
	int ID;
	float X;
	float Y;
	vector<int>neighbours;

};
struct Mesh {
	int num;
	float m_xub;
	float m_xlb;
	float m_yub;
	float m_ylb;
	float m_xlbc;
	float m_xubc;
	float m_yubc;
	float m_ylbc;
	float m_y;
	float s_rxlb;
	float s_rylb;
	float s_rxub;
	float s_ryub;
	int X1;
	int Y1;
	vector <int> memb;


};
vector<Particle>p;
vector<Mesh>m;
void createMesh(int& i, int& j, int& c, float& l_mb, float& l_m, float& Dx, int& num_grid)
{
	for (i = 0; i < num_grid; i++)
	{
		for (j = 0; j < num_grid; j++)
		{
			Mesh mm;
			mm.num = c; //ID of cell
			mm.m_xlb = (i * l_m - l_mb) * Dx; //uper x coordinate of each cell
			mm.m_xub = (i * l_m - l_mb) * Dx + l_m * Dx; // lower x coordinate
			mm.m_ylb = (j * l_m - l_mb) * Dx;// upper y coordinate of each cell
			mm.m_yub = (j * l_m - l_mb) * Dx + l_m * Dx; // lower y coordinate
			mm.m_xlbc = i;
			mm.m_xubc = i + 1;
			mm.m_ylbc = j;
			mm.m_yubc = j + 1;
			mm.s_rxlb = ((i - 1) * l_m - l_mb) * Dx;
			mm.s_rylb = ((j - 1) * l_m - l_mb) * Dx;
			mm.s_rxub = ((i + 1) * l_m - l_mb) * Dx + l_m * Dx;
			mm.s_ryub = ((j + 1) * l_m - l_mb) * Dx + l_m * Dx;
			c = c + 1;
			mm.X1 = i;
			mm.Y1 = j;

			m.push_back(mm);
		}
	}

}


void populatePartStruct(int& id, float& L, float& N, int& num_grid, float& Dx, int& v, int& a, int& b, int& c)
{
	Particle pp;
	for (int a = 0; a < N; a++)
	{
		for (int b = 0; b < N; b++)
		{

			pp.ID = id;
			id++;
			pp.X = a * Dx;
			pp.Y = b * Dx;

			p.push_back(pp);
		}
	}
	for (int h = 0; h < m.size(); h++)
	{
		for (int v = 0; v < p.size(); v++)
		{

			if (m[h].m_xlb < p[v].X && p[v].X < m[h].m_xub && m[h].m_ylb < p[v].Y && p[v].Y < m[h].m_yub)
			{
				m[h].memb.push_back(p[v].ID);
			}

		}
	}
}

void Neighboursearch(int w, int u, int t, int t_1, int t_2, int z, int z_1, int z_2, float distance1, float distance2, float distance3, int num_grid, float Dx)
{
	for (w = 0; w < m.size(); w++)
	{
		for (u = 0; u < m[w].memb.size(); u++)
		{
			if (int(w - (num_grid - 1)) % num_grid == 0)
			{


				for (int t = (m[w].num - (num_grid + 1)); t < (m[w].num - num_grid); t++)
				{
					if (t < 0 || t >(num_grid * num_grid) - 1)
					{
						continue;
					}

					for (int z = 0; z < m[t].memb.size(); z++)
					{
						float distance1 = sqrt((p[m[w].memb[u]].X - p[m[t].memb[z]].X) * (p[m[w].memb[u]].X - p[m[t].memb[z]].X) + (p[m[w].memb[u]].Y - p[m[t].memb[z]].Y) * (p[m[w].memb[u]].Y - p[m[t].memb[z]].Y));
						if (distance1 < 2 * Dx && distance1>0)
						{
							p[m[w].memb[u]].neighbours.push_back(p[m[t].memb[z]].ID);
						}
					}

				}
			}

			else if (w % (num_grid-1) == 0)
			{
				for (int t = (m[w].num - (num_grid)); t < (m[w].num - (num_grid - 1)); t++)
				{

					if (t < 0 || t >(num_grid * num_grid) - 1)
					{
						continue;


						for (int z = 0; z < m[t].memb.size(); z++)
						{
							distance1 = sqrt((p[m[w].memb[u]].X - p[m[t].memb[z]].X) * (p[m[w].memb[u]].X - p[m[t].memb[z]].X) + (p[m[w].memb[u]].Y - p[m[t].memb[z]].Y) * (p[m[w].memb[u]].Y - p[m[t].memb[z]].Y));
							if (distance1 < 2 * Dx && distance1>0)
							{
								p[m[w].memb[u]].neighbours.push_back(p[m[t].memb[z]].ID);
							}

						}
					}
				}

			}
			else
			{
				for (int t = (m[w].num - (num_grid + 1)); t < (m[w].num - (num_grid - 1)); t++)
				{

					if (t < 0 || t>(num_grid * num_grid) - 1)
					{
						continue;
					}

					for (int z = 0; z < m[t].memb.size(); z++)
					{
						distance1 = sqrt((p[m[w].memb[u]].X - p[m[t].memb[z]].X) * (p[m[w].memb[u]].X - p[m[t].memb[z]].X) + (p[m[w].memb[u]].Y - p[m[t].memb[z]].Y) * (p[m[w].memb[u]].Y - p[m[t].memb[z]].Y));
						if (distance1 < 2 * Dx && distance1>0)
						{
							p[m[w].memb[u]].neighbours.push_back(p[m[t].memb[z]].ID);
						}
					}

				}

			}


			if (w % (num_grid-1) == 0)
			{
				
				for (int t_1 = (m[w].num); t_1 < (m[w].num + 1); t_1++)
				{
					
					

					for (int z_1 = 0; z_1 < m[t_1].memb.size(); z_1++)
					{
						

						distance2 = sqrt((p[m[w].memb[u]].X - p[m[t_1].memb[z_1]].X) * (p[m[w].memb[u]].X - p[m[t_1].memb[z_1]].X) + (p[m[w].memb[u]].Y - p[m[t_1].memb[z_1]].Y) * (p[m[w].memb[u]].Y - p[m[t_1].memb[z_1]].Y));
						if (distance2 < 2 * Dx && distance2>0)
						{
							
							p[m[w].memb[u]].neighbours.push_back(p[m[t_1].memb[z_1]].ID);

							
						}
					}
				}
			}
			else if ((w - (num_grid - 1)) % num_grid == 0)
			{
				for (int t_1 = (m[w].num - 1); t_1 < (m[w].num); t_1++)
				{

					for (int z_1 = 0; z_1 < m[t_1].memb.size(); z_1++)
					{

						distance2 = sqrt((p[m[w].memb[u]].X - p[m[t_1].memb[z_1]].X) * (p[m[w].memb[u]].X - p[m[t_1].memb[z_1]].X) + (p[m[w].memb[u]].Y - p[m[t_1].memb[z_1]].Y) * (p[m[w].memb[u]].Y - p[m[t_1].memb[z_1]].Y));
						if (distance2 < 2 * Dx && distance2>0)
						{
							p[m[w].memb[u]].neighbours.push_back(p[m[t_1].memb[z_1]].ID);
						}
					}
				}
			}
			else
			{
				for (int t_1 = (m[w].num - 1); t_1 < (m[w].num + 1); t_1++)
				{
					if (t_1 < 0 || t_1 >(num_grid * num_grid) - 1)
					{
						continue;
					}

					for (int z_1 = 0; z_1 < m[t_1].memb.size(); z_1++)
					{

						distance2 = sqrt((p[m[w].memb[u]].X - p[m[t_1].memb[z_1]].X) * (p[m[w].memb[u]].X - p[m[t_1].memb[z_1]].X) + (p[m[w].memb[u]].Y - p[m[t_1].memb[z_1]].Y) * (p[m[w].memb[u]].Y - p[m[t_1].memb[z_1]].Y));
						if (distance2 < 2 * Dx && distance2>0)
						{
							p[m[w].memb[u]].neighbours.push_back(p[m[t_1].memb[z_1]].ID);
						}
					}

				}
			}

			if (int(w - (num_grid - 1) % num_grid == 0))
			{

				for (int t_2 = (m[w].num + (num_grid - 1)); t_2 < (m[w].num + (num_grid)); t_2++)
				{


					for (int z_2 = 0; z_2 < m[t_2].memb.size(); z_2++)
					{
						distance3 = sqrt((p[m[w].memb[u]].X - p[m[t_2].memb[z_2]].X) * (p[m[w].memb[z_2]].X - p[m[t_2].memb[z_2]].X) + (p[m[w].memb[u]].Y - p[m[t_2].memb[z_2]].Y) * (p[m[w].memb[u]].Y - p[m[t_2].memb[z_2]].Y));
						if (distance3 < 2 * Dx && distance3>0)
						{
							p[m[w].memb[u]].neighbours.push_back(p[m[t_2].memb[z_2]].ID);
						}
					}

				}
			}
			else if (w % num_grid == 0)
			{
				for (int t_2 = (m[w].num + (num_grid)); t_2 < (m[w].num + (num_grid + 1)); t_2++)
				{

					if (t_2 < 0 || t_2 >(num_grid * num_grid - 1))
					{
						continue;
					}

					for (int z_2 = 0; z_2 < m[t_2].memb.size(); z_2++)
					{
						distance3 = sqrt((p[m[w].memb[u]].X - p[m[t_2].memb[z_2]].X) * (p[m[w].memb[z_2]].X - p[m[t_2].memb[z_2]].X) + (p[m[w].memb[u]].Y - p[m[t_2].memb[z_2]].Y) * (p[m[w].memb[u]].Y - p[m[t_2].memb[z_2]].Y));
						if (distance3 < 2 * Dx && distance3 > 0)
						{
							p[m[w].memb[u]].neighbours.push_back(p[m[t_2].memb[z_2]].ID);
						}
					}

				}

			}
			else
			{
				for (int t_2 = (m[w].num + (num_grid - 1)); t_2 < (m[w].num + (num_grid + 1)); t_2++)
				{

					if (t_2 < 0 || t_2 >(num_grid * num_grid - 1))
					{
						continue;
					}

					for (int z_2 = 0; z_2 < m[t_2].memb.size(); z_2++)
					{
						distance3 = sqrt((p[m[w].memb[u]].X - p[m[t_2].memb[z_2]].X) * (p[m[w].memb[z_2]].X - p[m[t_2].memb[z_2]].X) + (p[m[w].memb[u]].Y - p[m[t_2].memb[z_2]].Y) * (p[m[w].memb[u]].Y - p[m[t_2].memb[z_2]].Y));
						if (distance3 < 2 * Dx && distance3>0)
						{
							p[m[w].memb[u]].neighbours.push_back(p[m[t_2].memb[z_2]].ID);
						}
					}

				}

			}
		}
	}
}

int main()
{
	int id = 0;
	float L = 1;
	float N = 100;
	int num_grid = 20;

	int i, j, a, b, ID, v, h;
	int c = 1;
	float x = 0;
	float y = 0;
	float Dx = L / (N - 1); //define the distance btw two particles (here the x and y are the same)
	float l_m = 5;//number of coordinates in x and y direction
	float l_mb = 0.5; //buffer on each side of mesh
	float dim_x = l_m * Dx;
	float dim_y = l_m * Dx;
	float x_m = 0;
	float y_m = 0;
	int w = 0;
	int u = 0;
	int t = 0;
	int t_1 = 0;
	int t_2 = 0;
	int z = 0;
	int z_1 = 0;
	int z_2 = 0;
	float distance1 = 0;
	float distance2 = 0;
	float distance3 = 0;
	int q;

	createMesh(i, j, c, l_mb, l_m, Dx, num_grid);
	populatePartStruct(id, L, N, num_grid, Dx, v, a, b, c);
	Neighboursearch(w, u, t, t_1, t_2, z, z_1, z_2, distance1, distance2, distance3, num_grid, Dx);



	for (int d = 0; d < m.size(); d++)
	{
		std::cout << "ID: " <<m[d].num << endl;
		
		for (int e = 0; e < m[d].memb.size(); e++)
		{
			std::cout << m[d].memb[e] << " ";
		}
		cout << endl;
	}
}

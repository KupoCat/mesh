#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>
#include <string>
#include <iostream>
#include <map>
#include <igl/copyleft/cgal/wire_mesh.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <imgui/imgui.h>
#include <windows.h>

using namespace std;

char c;

template<typename T>
void printvec(std::vector<T> vec)
{
	for (auto t : vec)
	{
		cout << t << ' ';
	}
	cout << endl;
}

vector<int> range(int end, int start = 0, int step = 1)
{
	vector<int> out;
	for (int i = start; i < end; i += step)
	{
		out.push_back(i);
	}
	return out;
}

void read_cubic_mesh(string path, vector<vector<double>>& vertices, vector<vector<int>>& cubes)
{
	enum mode { ver = 1, cube };
	std::ifstream inFile;
	inFile.open(path);
	int mode = 0;
	int v, c;
	for (string line; std::getline(inFile, line); )
	{
		if (line.compare("Vertices") == 0)
		{
			std::getline(inFile, line);
			v = std::stoi(line);
			cout << "V: " << v << endl;
			mode = ver;
			continue;
		}
		if (line.compare("Hexahedra") == 0)
		{
			std::getline(inFile, line);
			c = std::stoi(line);
			cout << "C: " << c << endl;
			mode = cube;
			continue;
		}
		if (line.compare("End") == 0)
		{
			cout << "Finished reading from file. V: " << v << " C: " << c << endl;
			return;
		}
		if (mode == ver)
		{
			std::stringstream iss(line);
			double poss;
			std::vector<double> ver;
			while (iss >> poss)
				ver.push_back(poss);
			ver.pop_back();
			vertices.push_back(ver);
			v--;
		}
		if (mode == cube)
		{
			std::stringstream iss(line);
			int ver;
			std::vector<int> cube;
			while (iss >> ver)
				cube.push_back(ver - 1);
			cube.pop_back();
			cubes.push_back(cube);
			c--;
		}
	}
}

void get_edges_from_cube_mesh(vector<vector<int>>& cubes, vector<vector<int>>& edges)
{
	for (auto cube : cubes)
	{
		for (auto i : range(3))
		{
			vector<int> edge;
			edge.push_back(cube[i]);
			edge.push_back(cube[i + 1]);
			edges.push_back(edge);
			edge.clear();
			edge.push_back(cube[i]);
			edge.push_back(cube[i + 4]);
			edges.push_back(edge);
			edge.clear();
			edge.push_back(cube[i + 4]);
			edge.push_back(cube[i + 5]);
			edges.push_back(edge);
		}
		vector<int> edge;
		edge.push_back(cube[3]);
		edge.push_back(cube[0]);
		edges.push_back(edge);
		edge.clear();
		edge.push_back(cube[3]);
		edge.push_back(cube[7]);
		edges.push_back(edge);
		edge.clear();
		edge.push_back(cube[7]);
		edge.push_back(cube[4]);
		edges.push_back(edge);
	}
}

void make_avg_ver(int a, int b, vector<vector<double>>& vertices)
{
	vector<double> ver;
	for (int i : range(3))
	{
		double tmp = (vertices[a][i] * 0.8) + (vertices[b][i] * 0.2);
		ver.push_back(tmp);
	}
	vertices.push_back(ver);
	ver.clear();
	for (int i : range(3))
	{
		double tmp = (vertices[b][i] * 0.8) + (vertices[a][i] * 0.2);
		ver.push_back(tmp);
	}
	vertices.push_back(ver);
}

void make_cell(vector<vector<double>>& vertices, vector<vector<int>>& cubes)
{
	// gets cubes and vertices and generates smaller cubes inside the bigger cubes //
	vector<vector<int>> tmp = cubes;
	for (auto cube : cubes)
	{
		vector<int> new_cube;
		int ind = vertices.size() - 1;

		make_avg_ver(cube[0], cube[6], vertices);
		make_avg_ver(cube[1], cube[7], vertices);
		make_avg_ver(cube[2], cube[4], vertices);
		make_avg_ver(cube[3], cube[5], vertices);

		for (int i : range(4))
		{
			new_cube.push_back(ind + (2 * i) + 1);
		}
		new_cube.push_back(ind + 6);
		new_cube.push_back(ind + 8);
		new_cube.push_back(ind + 2);
		new_cube.push_back(ind + 4);
		tmp.push_back(new_cube);
	}
	cubes = tmp;
}

void connect_cells(vector<vector<int>>& edges)
{
	int edges_size = edges.size() / 24; // number of original cubes = 12(number of edges in cube) * 2(half of the cubes are secondary)
	vector<vector<int>> tmp = edges;

	for (int i : range(edges_size)) // loop over cubes
	{
		for (int j : range(4)) // loop over sides
		{
			vector<int> edge;
			//top connector
			int a = edges[(i * 12) + ((j * 3) + 1)][0]; // cube i, side j, middle edge (|), starting point
			int b = edges[(i * 12) + ((edges_size * 12)) + ((j * 3) + 1)][0];// cube i of the secondary cubes, side j, middle edge (|), starting point
			edge.push_back(a);
			edge.push_back(b);
			tmp.push_back(edge);
			edge.clear();
			//bottom connector
			a = edges[(i * 12) + ((j * 3) + 1)][1];
			b = edges[(i * 12) + ((edges_size * 12)) + ((j * 3) + 1)][1];
			edge.push_back(a);
			edge.push_back(b);
			tmp.push_back(edge);
		}
	}
	set<vector<int>> s(tmp.begin(), tmp.end());
	edges = vector<vector<int>>(s.begin(), s.end());
}

double interpolate(double a, double b, double frac) // a, b - two points on a line. frac - fraction at which to interpolate
{
	return (a * frac) + (b * (1 - frac));
}

vector<vector<double>> convert_cube(vector<int>& source, vector<vector<double>>& from_ver, vector<int> target, vector<vector<double>> to_ver) // returns a vector of vertices of the interpolated cube
{
	vector<vector<double>> small;
	vector<vector<double>> big;
	vector<vector<double>> out;
	for (auto ind : source)
	{
		small.push_back(from_ver[ind]);
	}
	for (auto ind : target)
	{
		big.push_back(to_ver[ind]);
	}
	
	int i = 0;
	for (auto ver : small)
	{
		vector<double> new_ver;
		for (int dim : range(3)) // interpolate for each dim
		{
			//first interpolation
			double a = interpolate(big[0][dim], big[1][dim], ver[0]);
			double b = interpolate(big[3][dim], big[2][dim], ver[0]);
			double c = interpolate(big[4][dim], big[5][dim], ver[0]);
			double d = interpolate(big[7][dim], big[6][dim], ver[0]);
			//second interpolation
			double a_ = interpolate(a, b, ver[1]);
			double b_ = interpolate(c, d, ver[1]);
			//third interpolation
			double x = interpolate(a_, b_, ver[2]);
			new_ver.push_back(x);
		}
		out.push_back(new_ver);
	}
	return out;
}

void import_cell(vector<vector<double>>& M_vertices, vector<vector<int>>& M_cubes, vector<vector<double>>& C_vertices, vector<vector<int>>& C_cubes)
{
	vector<vector<int>> tmp_cubes;
	vector<vector<double>> tmp_vers;
	int jump = M_vertices.size();
	for (auto M_cube : M_cubes)
	{
		for (auto C_cube : C_cubes)
		{
			vector<vector<double>> new_vers = convert_cube(C_cube, C_vertices, M_cube, M_vertices);
			vector<int> new_cube = range(jump + 8, jump);
			jump += 8;
			tmp_cubes.push_back(new_cube);
			tmp_vers.insert(tmp_vers.end(), new_vers.begin(), new_vers.end());
		}
	}

	M_vertices.insert(M_vertices.end(), tmp_vers.begin(), tmp_vers.end());
	M_cubes.insert(M_cubes.end(), tmp_cubes.begin(), tmp_cubes.end());
}

string ExePath() {
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	string::size_type pos = string(buffer).find_last_of("\\/");
	return string(buffer).substr(0, pos);
}

vector<vector<int>> mesh_to_wire(vector<vector<double>> vertices, vector<vector<int>> cubes, string full_path, float thickness, int faces, bool clean)
{
	vector<vector<int>> edges;

	get_edges_from_cube_mesh(cubes, edges);

	cout << "got edges" << endl;
	// in
	Eigen::Matrix<double, Eigen::Dynamic, 3> MatVertices;
	Eigen::Matrix<int, Eigen::Dynamic, 2> MatEdges;

	igl::list_to_matrix(vertices, MatVertices);
	igl::list_to_matrix(edges, MatEdges);

	// out
	Eigen::Matrix<int, Eigen::Dynamic, 3> outFaces;
	Eigen::Matrix<double, Eigen::Dynamic, 3> outVertecies;
	Eigen::Matrix<int, Eigen::Dynamic, 1> J;

	igl::copyleft::cgal::wire_mesh(MatVertices, MatEdges, thickness, faces, clean, outVertecies, outFaces, J);
	igl::writeOBJ(full_path, outVertecies, outFaces); // write wiremesh to obj

	return edges;
}

int main(int argc, char *argv[])
{
	std::string base_path = ExePath() + std::string("\\..\\..\\Data\\");

	vector<vector<double>> Model_ver;
	vector<vector<int>> Model_cubes;
	vector<vector<double>> Cell_ver;
	vector<vector<int>> Cell_cubes;
	vector<vector<int>> ModelEdges;
	vector<vector<int>> CellEdges;

	float f = 0.01f;
	int i = 3;
	bool b = false;

	read_cubic_mesh(base_path + "bunny_1.mesh", Model_ver, Model_cubes);
	read_cubic_mesh(base_path + "cell.mesh", Cell_ver, Cell_cubes);

	//make_cell(inVertices, Cubes);

	import_cell(Cell_ver, Cell_cubes, Cell_ver, Cell_cubes);

	import_cell(Model_ver, Model_cubes, Cell_ver, Cell_cubes);

	// -------------------MODEL------------------- //

	ModelEdges = mesh_to_wire(Model_ver, Model_cubes, base_path + "\\model.obj", f, i, b);
	
	// -------------------CELL------------------- //

	CellEdges = mesh_to_wire(Cell_ver, Cell_cubes, base_path + "\\cell.obj", f, i, b);

	// --------------------MENU-------------------- //
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	viewer.load_mesh_from_file(base_path + "\\cell.obj");
	viewer.load_mesh_from_file(base_path + "\\model.obj");

	int left_view, right_view;
	int cell_id = viewer.data_list[0].id, model_id = viewer.data_list[1].id;

	cout << "split viewer" << endl;


	viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
	{
		viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
		left_view = viewer.core_list[0].id;
		right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));
		cout << viewer.core_list.size() << endl;
		return false;
	};

	viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
	{
		if (key == GLFW_KEY_SPACE)
		{
			// by default, when a core is appended, all loaded meshes will be displayed in that core
			// displaying can be controlled by changing viewer.coreDataPairs
			viewer.data(cell_id).set_visible(false, left_view);
			viewer.data(model_id).set_visible(false, right_view);
		}
		return false;
	};

	viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h) {
		v.core(left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
		v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h);
		return true;
	};

	

	vector<double> ver;
	ver.push_back(0);
	ver.push_back(0);
	ver.push_back(0);

	menu.callback_draw_viewer_menu = [&]()
	{
		ImGui::InputFloat("Thickness", &f);
		ImGui::InputInt("Polygon size", &i);
		ImGui::Checkbox("Clean Mesh", &b);

		if (ImGui::Button("Reload wire mesh"))
		{
			mesh_to_wire(Model_ver, Model_cubes, base_path + "\\model.obj", f, i, b);
			viewer.data(model_id).clear();

			Eigen::Matrix<int, Eigen::Dynamic, 3> F;
			Eigen::Matrix<double, Eigen::Dynamic, 3> V;

			igl::readOBJ(base_path + "\\model.obj", V, F);

			viewer.data(model_id).set_mesh(V, F);
		}
		ImGui::NewLine();

		ImGui::Text("Edit Cell");
		ImGui::Text("---------------------");
		ImGui::InputDouble("X: ", &ver[0]);
		ImGui::InputDouble("Y: ", &ver[1]);
		ImGui::InputDouble("Z: ", &ver[2]);
		if (ImGui::Button("Add Vertex"))
		{

		}
	};
	
	
	cout << "done" << endl;
	viewer.launch();
}
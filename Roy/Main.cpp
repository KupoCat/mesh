#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>
#include <string>
#include <iostream>
#include <map>
#include <igl/copyleft/cgal/wire_mesh.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <imgui/imgui.h>

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

int main(int argc, char *argv[])
{

	vector<vector<double>> inVertices;
	vector<vector<int>> Cubes;
	vector<vector<int>> inEdges;

	read_cubic_mesh("C:\\Users\\tomer\\Desktop\\python\\Mesh to Obj\\hand.mesh", inVertices, Cubes);

	//Cubes = vector<vector<int>>(Cubes.begin(), Cubes.begin() + 10);

	make_cell(inVertices, Cubes);

	cout << "finished making cells" << endl;

	get_edges_from_cube_mesh(Cubes, inEdges);

	cout << "got edges" << endl;

	connect_cells(inEdges);

	cout << "connected cells" << endl;

	// in
	Eigen::Matrix<double, Eigen::Dynamic, 3> MatVertices;
	Eigen::Matrix<int, Eigen::Dynamic, 2> MatEdges;

	igl::list_to_matrix(inVertices, MatVertices);
	igl::list_to_matrix(inEdges, MatEdges);

	// out
	Eigen::Matrix<int, Eigen::Dynamic, 3> outFaces;
	Eigen::Matrix<double, Eigen::Dynamic, 3> outVertecies;
	Eigen::Matrix<int, Eigen::Dynamic, 1> J;

	Eigen::Matrix<int, Eigen::Dynamic, 3> cellFace;
	Eigen::Matrix<double, Eigen::Dynamic, 3> cellVer;

	igl::copyleft::cgal::wire_mesh(MatVertices, MatEdges, 0.01, 3, false, outVertecies, outFaces, J);
	cout << "Done" << endl;
	// --------------------MENU-------------------- //
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);
	std::string base_path = std::string("C:\\Users\\Tomer\\Desktop\\libigl\\tutorial\\data");
	//igl::writeOBJ(base_path + "\\hand.obj", outVertecies, outFaces);
	
	viewer.load_mesh_from_file(base_path + "\\cube.obj");
	viewer.load_mesh_from_file(base_path + "\\hand.obj");
	int left_view, right_view;
	int cell_id = viewer.data_list[0].id, model_id = viewer.data_list[1].id;


	//cout << viewer.data() << endl;

	viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
	{
		viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
		left_view = viewer.core_list[0].id;
		right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));
		return true;
	};
	cout << "split viewer" << endl;
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

	float f = 0.01f;
	int i = 3;
	bool b = false;
	cout << "loading menu" << endl;
	/*
	menu.callback_draw_viewer_menu = [&]()
	{
		ImGui::InputFloat("Thickness", &f);
		ImGui::InputInt("Polygon size", &i);
		ImGui::Checkbox("Clean Mesh", &b);

		if (ImGui::Button("Reload wire mesh"))
		{
			igl::copyleft::cgal::wire_mesh(MatVertices, MatEdges, f, i, b, outVertecies, outFaces, J); // wire mesh
			viewer.data().clear();
			viewer.data().set_mesh(outVertecies, outFaces); // display
		}
	};
	*/
	cout << "done" << endl;
	viewer.launch();
	cin >> c;
}
/*
int main(int argc, char * argv[])
{
	igl::opengl::glfw::Viewer viewer;

	viewer.load_mesh_from_file(std::string("C:\\Users\\Tomer\\Desktop\\libigl\\tutorial\\data") + "/cube.obj");
	viewer.load_mesh_from_file(std::string("C:\\Users\\Tomer\\Desktop\\libigl\\tutorial\\data") + "/sphere.obj");

	int left_view, right_view;
	int cube_id = viewer.data_list[0].id, sphere_id = viewer.data_list[1].id;
	viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
	{
		viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
		left_view = viewer.core_list[0].id;
		right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));
		return true;
	};

	viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
	{
		if (key == GLFW_KEY_SPACE)
		{
			// by default, when a core is appended, all loaded meshes will be displayed in that core
			// displaying can be controlled by changing viewer.coreDataPairs
			viewer.data(cube_id).set_visible(false, left_view);
			viewer.data(sphere_id).set_visible(false, right_view);
		}
		return false;
	};

	viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h) {
		v.core(left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
		v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h);
		return true;
	};

	viewer.launch();
	return EXIT_SUCCESS;
}
*/
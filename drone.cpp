// Project Identifier: 1761414855B69983BD8035097EFBD312EB0527F0

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>
#include <iomanip>
using namespace std;

class Coord {
    public: 
        int x = 0;
        int y = 0;
        string area;
        bool visited = false;
        double weight = std::numeric_limits<double>::infinity();
        size_t parent = 0;
};

class Drone {
    public:
        void get_options(int argc, char** argv);
        void read_input();
        void mst();
        void arbitrary_insertion();
        void fasttsp();
        void opttsp();
        void promising(size_t & start_ind, double & sum);
        void print();
        void two_opt(vector<size_t> & v);
        void genPerms(size_t permLength);
        bool promising(size_t permLength);
        string mode;
        size_t v = 0;

    private:
        vector<Coord> vertices;
        vector<size_t> path;
        double curCost = 0;
        vector<size_t> best_path;
        double best_cost = 0;
        double sum = 0;
};

double mst_distance(Coord c1, Coord c2) {
    if ((c1.area != c2.area) && (c1.area != "border" && c2.area != "border")) {
        return std::numeric_limits<double>::infinity();
    }
    int x_ = c2.x - c1.x;
    int y_= c2.y - c1.y;
    double x = static_cast<double>(x_);
    double y = static_cast<double>(y_);
    return (x * x) + (y * y);
}

double fasttsp_distance(Coord c1, Coord c2) {
    int x_ = c2.x - c1.x;
    int y_ = c2.y - c1.y;
    double x = static_cast<double>(x_);
    double y = static_cast<double>(y_);
    return sqrt((x * x) + (y * y));
}

double distance(Coord c1, Coord c2) {
    int x_ = c2.x - c1.x;
    int y_ = c2.y - c1.y;
    double x = static_cast<double>(x_);
    double y = static_cast<double>(y_);
    return (x * x) + (y * y);
}

int main(int argc, char** argv) {
    std::ios_base::sync_with_stdio(false);
    Drone d;
    cout << fixed << showpoint << setprecision(2) << boolalpha;

    d.get_options(argc, argv);
    d.read_input();

    if (d.mode == "MST") {
        d.mst();
    }else if (d.mode == "FASTTSP") {
        d.arbitrary_insertion();
    }else if (d.mode == "OPTTSP") {
        d.opttsp();
    }
    d.print();

    return 0;
}

void Drone::read_input() {
    cin >> v;
    vertices.reserve(v);

    for (size_t i = 0; i < v; ++i) {
        Coord c;
        int x_ = 0;
        int y_ = 0;
        cin >> x_ >> y_;
        c.x = x_;
        c.y = y_;
        if (x_ < 0 && y_ < 0) {
            c.area = "medical";
        }else if (x_ == 0 || y_ == 0) {
            c.area = "border";
        }else {
            c.area = "main";
        }
        vertices.push_back(c);
    }
}


void Drone::mst() {
    size_t min_ind = 0;
    vertices[0].weight = 0;
    for (size_t count = 0; count < v; ++count) {
        double min = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < vertices.size(); ++i) {
            if (vertices[i].weight < min && vertices[i].visited == false) {
                min_ind = i;
                min = vertices[i].weight;
            }
        }
        vertices[min_ind].visited = true;
        sum += vertices[min_ind].weight;
        for (size_t j = 0; j < vertices.size(); ++j) {
            if (vertices[j].visited == false) {
                double d = mst_distance(vertices[min_ind], vertices[j]);
                if (d < (vertices[j].weight * vertices[j].weight)) {
                    vertices[j].weight = sqrt(d);
                    vertices[j].parent = min_ind;
                }
            }
        }
    }
}

void reverse(vector<size_t> & v, int begin, int end) {
    while (end - begin > 0) {
        size_t e = size_t(end);
        size_t b = size_t(begin);
        size_t temp = v[b % v.size()];
        v[b % v.size()] = v[e % v.size()];
		v[e % v.size()] = temp;
		begin++;
		end--;
    }
}

void Drone::two_opt(vector<size_t> & v) {
    for (size_t i = 1; i < v.size() - 2; ++i) {
        for (size_t j = i + 2; j < v.size(); ++j) {
                size_t next = (j + 1) % v.size();
                double orig_dist = fasttsp_distance(vertices[v[i]], vertices[v[(i+1)]]) + fasttsp_distance(vertices[v[j]], vertices[v[next]]);
                double new_dist = fasttsp_distance(vertices[v[i]], vertices[v[j]]) + fasttsp_distance(vertices[v[(i+1)]], vertices[v[(next)]]);
                if (new_dist < orig_dist) {
                    std::reverse(v.begin() + (int)i + 1, v.begin() + (int)j + 1);
                    sum = sum - orig_dist + new_dist;
                }
        }
    }
}

void Drone::arbitrary_insertion() {
    path.reserve(v + 1);
    vertices[0].visited = true;
    vertices[1].visited = true;
    path.push_back(0);
    path.push_back(1);
    path.push_back(0);

    sum += 2 * fasttsp_distance(vertices[path[1]], vertices[0]);

    for (size_t i = 2; i < vertices.size(); ++i) {
        double min_dist = std::numeric_limits<double>::infinity();
        size_t k = i;
        size_t min_ind = 0;
        for (size_t j = 0; j < path.size(); ++j) {
                size_t next = j + 1;
                if (next == path.size()) {
                    next = 0;
                }
                double dist1 = fasttsp_distance(vertices[path[j]], vertices[k]);
                double dist2 = fasttsp_distance(vertices[k], vertices[path[next]]);
                double new_dist = dist1 + dist2 - fasttsp_distance(vertices[path[j]], vertices[path[next]]);;
                if (new_dist < min_dist) {
                    min_dist = new_dist;
                    min_ind = next;
                }
        }
        sum += min_dist;
        path.insert(path.begin() + (int)min_ind, k);
    }
    path.pop_back();
    sum += fasttsp_distance(vertices[path[path.size() - 1]], vertices[0]);
}

void Drone::fasttsp() {
    path.reserve(v);
    size_t curr_ind = 0;
    vertices[0].visited = true;
    path.push_back(curr_ind);

    while (path.size() < v) {
        double min_dist = std::numeric_limits<double>::infinity();
        size_t min_ind = curr_ind;
        for (size_t i = 0; i < v; ++i) {
            if (vertices[i].visited == false) {
                double d = distance(vertices[curr_ind], vertices[i]);
                if (d < (min_dist * min_dist)) {
                    min_ind = i;
                    min_dist = sqrt(d);
                }
            }
        }
        vertices[min_ind].visited = true;
        curr_ind = min_ind;
        path.push_back(curr_ind);
        sum += min_dist;
    }
    sum += fasttsp_distance(vertices[curr_ind], vertices[0]);

    two_opt(path);
}


void Drone::print() {
    if (mode == "MST") {
        cout << sum << "\n";
        for (size_t i = 1; i < vertices.size(); ++i) {
            if (i < vertices[i].parent) {
                cout << i << " " << vertices[i].parent << "\n";
            }else {
                cout << vertices[i].parent << " " << i << "\n";
            }
        }
    }else if (mode == "FASTTSP") {
        cout << sum << "\n";
        for (size_t i = 0; i < path.size(); ++i) {
            cout << path[i] << " ";
        }
    }else if (mode == "OPTTSP") {
        cout << best_cost << "\n";
        for (size_t i = 0; i < best_path.size(); ++i) {
            cout << best_path[i] << " ";
        }
    }
}

bool Drone::promising(size_t permLength) {
    for (size_t i = permLength; i < path.size(); ++i) {
        vertices[path[i]].visited = false;
        vertices[path[i]].weight = std::numeric_limits<double>::infinity();
    }

    size_t min_ind = path[permLength];
    double mstCost = 0;
    vertices[min_ind].weight = 0;
    for (size_t count = permLength; count < best_path.size(); ++count) {
        double min = std::numeric_limits<double>::infinity();
        for (size_t i = permLength; i < best_path.size(); ++i) {
            if (vertices[path[i]].weight < min && vertices[path[i]].visited == false) {
                min_ind = path[i];
                min = vertices[path[i]].weight;
            }
        }
        vertices[min_ind].visited = true;
        mstCost += min;
        for (size_t j = permLength; j < best_path.size(); ++j) {
            if (vertices[path[j]].visited == false) {
                double d = fasttsp_distance(vertices[min_ind], vertices[path[j]]);
                if (d < vertices[path[j]].weight) {
                    vertices[path[j]].weight = d;
                }
            }
        }
    }

    double arm1Len = std::numeric_limits<double>::infinity();
    double arm2Len = std::numeric_limits<double>::infinity();
    size_t begin = 0;
    size_t end = path[permLength - 1];

    for (size_t i = permLength; i < path.size(); ++i) {
        double d = distance(vertices[begin], vertices[path[i]]);
        if (d < (arm1Len * arm1Len)) {
            arm1Len = sqrt(d);
        }
        d = distance(vertices[end], vertices[path[i]]);
        if (d < (arm2Len * arm2Len)) {
            arm2Len = sqrt(d);
        }
    }

    double totalEst = arm1Len + arm2Len + mstCost + curCost;
    bool promise = false;
    if (totalEst < best_cost) {
        promise = true;
    }

    return promise;
}

void Drone::genPerms(size_t permLength) {
    if (permLength == path.size()) {
        double d = fasttsp_distance(vertices[path[permLength - 1]], vertices[path[0]]);
        if (curCost + d < best_cost) {
            best_cost = curCost + d;
            best_path = path;
        }
        return;
    }  // if ..complete path

    if (!promising(permLength)) {
        return;
    }  // if ..not promising

    for (size_t i = permLength; i < path.size(); ++i) {
        swap(path[permLength], path[i]);
        curCost += fasttsp_distance(vertices[path[permLength - 1]], vertices[path[permLength]]);
        genPerms(permLength + 1);
        curCost -= fasttsp_distance(vertices[path[permLength - 1]], vertices[path[permLength]]);
        swap(path[permLength], path[i]);
    }  // for ..unpermuted elements
}  // genPerms()

void Drone::opttsp() {
    arbitrary_insertion();
    best_path = path;
    best_cost = sum;
    genPerms(1);
}

// Read and process command line options.
void Drone::get_options(int argc, char** argv) {
    int option_index = 0, option = 0;
    
    // Don't display getopt error messages about options
    opterr = false;

    /*

        TODO: Add the remaining elements into the longOpts array.

    */
    // use getopt to find command line options
    struct option longOpts[] = {{ "mode", required_argument, nullptr, 'm' },
                                { "help", no_argument, nullptr, 'h' },
                                { nullptr, 0, nullptr, '\0' }};
    
    /*

        TODO: Add the remaining chars to the option string in
                the while loop conditional (currently contains "p:h").
                Options with required_argument (print) need a colon after the
                char, options with no_argument do not (help).

    */
    while ((option = getopt_long(argc, argv, "m:h", longOpts, &option_index)) != -1) {
        switch (option) {
            case 'm': {
                string str(optarg);
                mode = str;
                break; 
            }    
    
            case 'h':
                std::cout << "This program reads a CSV file that contains song names,\n"
                          << "the artist who wrote them, and the number of plays each song\n"
                          <<  "has on Spotify.  It then outputs the number of songs specified\n"
                          <<  "in the command line arguments (the print option), which\n"
                          <<  "defaults to 2, sorted by the option specified (one of name,\n"
                          <<  "artist, or listens).\n"
                          <<  "Usage: \'./project0\n\t[--listens | -l]\n"
                          <<                      "\t[--name | -n]\n"
                          <<                      "\t[--artist | -a]\n"
                          <<                      "\t[--print | -p] <# of songs to print>\n"
                          <<                      "\t[--help | -h]\n"
                          <<                      "\t< <CSV Music File>\'" << std::endl;
                exit(0);
            
            default:
                cerr << "Unknown command line option" << endl;
                exit(0);
        }
    }
}

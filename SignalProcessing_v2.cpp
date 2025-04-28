#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <tuple>
#include <cstdlib>

using namespace std;

const int run_start = 1;
const int run_stop = 1;


string file_name = "C:\\code_root\\data_in\\f1";
const int events_per_run = 1000;
const int points_per_event = 9999;

typedef enum {
	CHANNEL_SLOW = 2, 
	CHANNEL_FAST = 5,
}CHANNEL;

typedef struct {
	int x;
	int y;
}Point;

struct Events{
	int event;
	int run;

	bool operator<(const Events& other) const {
        return tie(run, event) < tie(other.run, other.event);
    }

};

vector<int> Read_by_runs(CHANNEL ch, int run){
	vector <int> signal_run;
	signal_run.resize(events_per_run*points_per_event);

	FILE *f = NULL;
	
	ostringstream _run, _ch;
	_run << run;
	_ch  << ch;
	string channel =  "__ch_" + _ch.str() + ".dat";
	string filepath = file_name + "\\run_" + _run.str() + channel;

	f = fopen(filepath.c_str(), "rb");

	if (!f) {
		cout << "Error opening file: " << filepath << endl;
		exit(EXIT_FAILURE);
	}

	fread(signal_run.data(), 2, events_per_run*points_per_event, f);
	fclose(f); 

	return signal_run;
} 

vector<Point> Turn_origin_to_point(CHANNEL ch, vector<int> signal){
	int coeff;
	if (ch == CHANNEL_SLOW) coeff = 1;
	if (ch == CHANNEL_FAST) coeff = -1;
	
	vector<Point> signal_points;
	signal_points.reserve(signal.size()); 
	Point point;
	for (int points = 0; points < signal.size(); points++){
		point.x = points;
		point.y = coeff*signal[points];
		signal_points.push_back(point);
	} 
	return signal_points;
}

vector<vector<Point>> Parsing_data_on_events(vector<Point> signal){
	vector<vector<Point>> signal_events;
	for (int event = 0; event < events_per_run; event++)
	{
		vector <Point> temp(signal.begin() + event*points_per_event, signal.begin() + (event+1)*points_per_event);
		signal_events.push_back(std::move(temp));
	}

	return signal_events;
}

map<Events, vector<Point>> Mapping_data_with_runs_events(int run, vector<vector<Point>> signal){
	map<Events, vector<Point>> mapped_signal;
	for(int event = 0; event < signal.size(); event++){
		Events run_event;
		run_event.run = run;
		run_event.event = event;
		mapped_signal.insert({run_event, signal[event]}); 
	}
	return mapped_signal;
}



int Signal_processing2()
{	
	
	vector<vector<Point>> signal;
	signal = Parsing_data_on_events(Turn_origin_to_point(CHANNEL_SLOW, Read_by_runs(CHANNEL_SLOW, run_start)));
	auto ssignal = Mapping_data_with_runs_events(run_start, signal);
	for (const auto& pair : ssignal){
		cout << pair.first.run << " " << pair.first.event << " " << pair.second[9999].y << endl; 
	}
	return 0;
}
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <math.h>

#include <TFile.h>
#include <TTree.h>

using namespace std;
namespace fs = std::filesystem; // Just an alias

typedef int16_t data_t; // NOTE: the binary data must be encoded into this type!

const int RUN_START = 1;
const int RUN_STOP = 1;

const fs::path DATA_DIR = "./data/";
const int EVENTS_PER_RUN = 1000;
const int POINTS_PER_EVENT = 9999;
const int SEC_PER_POINT = 16;
const int PERIOD = 10000;

typedef enum {
    CHANNEL_SLOW,
    CHANNEL_FAST,
} CHANNEL_TYPE;

// TODO: If this is just a collection of constants we can simply make make like this:
//       ```
//       const double THRESHOLD_SLOW = 40;
//       const double THRESHOLD_FAST = 20;
//       ```
//       Or if we intend to use it as a type, then, we have to change the type of a
//       `Channel.threshold` variable to `THRESHOLD_TYPE`
typedef enum { 
    THRESH_SLOW = 40,
    THRESH_FAST = 20,
} THRESHOLD_TYPE; 

typedef struct {
    CHANNEL_TYPE t;
    int n;
    double threshold; // See comment above `THRESHOLD_TYPE` enum
} Channel;

typedef struct {
    int x;
    data_t y;
} Point;

typedef vector<Point> Event;

struct RunEvent_K {
    int run;
    int event;

    bool operator<(const RunEvent_K& other) const {
        if (run == other.run) return event < other.event;
        return run < other.run;
    }
};

typedef struct {
    int left;
    int right;
} Borders;


// two trees for each channel which fill with struct Event_data
// Event_data include struct Peak_data as vector of peaks
/*
struct Peak_data{
	double_t area;
	Point max;
	Event coord;
	Peak_data() : area(0) {}
};

struct Event_data{
	double_t total_area;
	double_t start_fluct;
	vector <Peak_data> peaks;
	//TODO add smth i forget
	Event_data() : total_area(0) {}
};

void creat_tree(){
	TFile *file = new TFile("event.root", "RECREATE");
	TTree *Chan_slow = new TTree("Chan_slow", "CHANNEL SLOW data");
	TTree *Chan_fast = new TTree("Chan_fast", "CHANNEL FAST data");
	// two trees of fastslow signals
	Event_data CHANNEL_SLOW;
	Event_data CHANNEL_FAST;
	Chan_slow->Branch("total_area", &CHANNEL_SLOW.total_area, "total_area/F");
	Chan_slow->Branch("peaks", &CHANNEL_SLOW.peaks); // This won't work directly; see below
    // Create branches for Channel 2
    Chan_fast->Branch("total_area", &CHANNEL_FAST.total_area, "total_area/F");
    Chan_fast->Branch("peaks", &CHANNEL_FAST.peaks);
	file->Write();
    file->Close();
}
*/

vector<data_t> read_data(fs::path filepath, size_t n_obj) {
    vector<data_t> data(n_obj);

    FILE* f = fopen(filepath.c_str(), "rb");
    if (!f) {
        cout << "Error opening file: " << filepath << endl;
        exit(EXIT_FAILURE);
    }
    fread(data.data(), sizeof(data_t), n_obj, f);
    fclose(f);

    return data;
}


// Read data from channel `ch` and run `run`
vector<data_t> read_data_by_run(Channel chan, int run) {
    char filename[32];
    sprintf(filename, "run_%d__ch_%d.dat", run, chan.n);
    fs::path filepath = DATA_DIR / filename;

    return read_data(filepath, POINTS_PER_EVENT * EVENTS_PER_RUN);
}


vector<Event> split_data_to_events(Channel chan, vector<data_t> data) {
    int coeff;
    if (chan.t == CHANNEL_SLOW) coeff = 1;
    else if (chan.t == CHANNEL_FAST) coeff = -1;

    vector<Event> events(EVENTS_PER_RUN);
    for (size_t i = 0; i < EVENTS_PER_RUN; i++) {
        Event event(POINTS_PER_EVENT);
        for (size_t j = 0; j < POINTS_PER_EVENT; j++) {
            event[j] = Point {
                .x = (int)j,
                .y = static_cast<data_t>(coeff*data[i * POINTS_PER_EVENT + j]),
            };
        }
        events[i] = event;
    }

    return events;
}

void find_peaks(Channel chan, Event data, Borders left_right = {}) {
    Borders peak_borders;
    Borders timezone_to_fp;
    double_t area, total_area;
    int temp_left, temp_right;
    int depth = 5;

    Event event = data;
    if (left_right.left != NULL ) {
        timezone_to_fp.left = 1;
        timezone_to_fp.right = data.size();
    }
    else {
        timezone_to_fp.left = left_right.left;
        timezone_to_fp.right = left_right.right;
    }
    if (timezone_to_fp.left < 0 || timezone_to_fp.right >event.size()) {
        cout << "out of limits" << endl;
        timezone_to_fp.left = 1;
        timezone_to_fp.right = event.size();
    }
    for (size_t i = timezone_to_fp.left + 1; i < timezone_to_fp.right; i++) {
        Point max = {0, 0};
        if(event[i-1].y < chan.threshold && event[i].y > chan.threshold) {
            peak_borders.left = i;
            while (event[i].y > chan.threshold/depth && peak_borders.left > 0) peak_borders.left--;
            peak_borders.right = i;
            while (event[i].y > chan.threshold/depth && peak_borders.right < POINTS_PER_EVENT) peak_borders.right++;

            if (peak_borders.left != temp_left && peak_borders.right != temp_right) {
                for (size_t i = peak_borders.left; i <= peak_borders.right; i++) {
                    area += SEC_PER_POINT*(event[i-1].y + event[i].y)/2;
                    if (max.y < event[i].y)
                        max = event[i];
                }
            }
            else area = 0;
            total_area += area;
        }
    }
    temp_left  = peak_borders.left;
    temp_right = peak_borders.right;

    //TODO obtain data to tree of peaks {area, max, peak_borders} mb its leaves, obtain data to tree of event {total area}
}

double_t check_start_fluct(int run, const vector<Event> data) {
    int startline = PERIOD/SEC_PER_POINT;
    double_t maxy_start = 0;
    for (size_t i = 0; i < startline; i++) {
        if ( maxy_start > data[run][i].y)
            maxy_start = data[run][i].y;
    }
    return maxy_start;
}

vector<Event> normalize_baseline(Channel chan, vector<Event> data) {
    double_t baseline_avr, baseline_sigma;
    int  sum, sum_sigma;
    int start_line = PERIOD/SEC_PER_POINT;
    double_t start_fluct;
    Borders timezone_to_fp = {.left = 0, .right = start_line};
    vector <Event> normal_signal = data;
    for (size_t i = 0; i < data.size(); i++) {
        sum = 0;
        sum_sigma = 0;
        // TODO obtain start_fluct to struct event
        start_fluct = check_start_fluct(i, normal_signal);
        for (size_t j = 0; j < start_line; j++) {
            sum += normal_signal[i][j].y;
        }
        baseline_avr = sum/start_line;

        for (size_t j = 0; j < POINTS_PER_EVENT; j++) {
            normal_signal[i][j].y -= baseline_avr;
        }
        for (size_t j = 0; j < start_line; j++) {
            sum_sigma += pow((normal_signal[i][j].y - baseline_avr), 2);
        }
        baseline_sigma = sqrt(sum_sigma/start_line);
        //TODO obtain data of events into tree of event: baseline_sigma, baseline_avr,
    }
    return normal_signal;
}

// TODO: Discuss intention of this function
// Get a map, where key consists of a pair of 2 values <run> <event>, and value is an event data
map<RunEvent_K, Event> get_mapped_events_for_run(int run, vector<Event> signal) {
    map<RunEvent_K, Event> mapped_signal;
    for (size_t i = 0; i < signal.size(); i++) {
        RunEvent_K run_event = {
            .run = run,
            .event = (int)i,
        };
        mapped_signal.insert({run_event, signal[i]});
    }
    return mapped_signal;
}

int SignalProcessing_v2() {
    Channel chan = {
        .t = CHANNEL_SLOW,
        .n = 2,
        .threshold = THRESH_SLOW,
        //for fast channel for ex it would be .n = 5 .treashold = 20;
    };
    vector<data_t> data = read_data_by_run(chan, RUN_START);
    vector<Event> events = split_data_to_events(chan, data);

    // TODO: Discuss the reason for getting a map of <run, event> -> event_data
    map<RunEvent_K, Event> mapped_events = get_mapped_events_for_run(RUN_START, events);
    for (auto pair : mapped_events) {
        cout << "[" << pair.first.run << ":" << pair.first.event << "]: " << pair.second[POINTS_PER_EVENT - 1].y << endl;
    }
    return 0;
}

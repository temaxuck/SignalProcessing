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

using namespace std;
namespace fs = std::filesystem; // Just an alias

const size_t DATA_OBJECT_SIZE = 2; // number of bytes per data object
typedef int16_t data_t; // sizeof(DATA_OBJECT_TYPE) == DATA_OBJECT_SIZE

const int RUN_START = 1;
const int RUN_STOP = 1;

const fs::path DATA_DIR = "./data/";
const int EVENTS_PER_RUN = 1000;
const int POINTS_PER_EVENT = 9999;

typedef enum {
    CHANNEL_SLOW,
    CHANNEL_FAST,
} CHANNEL_TYPE;

typedef struct {
    CHANNEL_TYPE t;
    int n;
} Channel;

typedef struct {
    int x;
    int y;
} Point;

typedef vector<Point> Event;

struct Events {
    int event;
    int run;

    bool operator<(const Events& other) const {
        return tie(run, event) < tie(other.run, other.event);
    }

};

vector<data_t> read_data(fs::path filepath, size_t n_obj) {
    vector<data_t> data(n_obj);

    FILE* f = fopen(filepath.c_str(), "rb");
    if (!f) {
        cout << "Error opening file: " << filepath << endl;
        exit(EXIT_FAILURE);
    }
    fread(data.data(), DATA_OBJECT_SIZE, n_obj, f);
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
                .y = coeff*data[i * POINTS_PER_EVENT + j],
            };
        }
        events[i] = event;
    }

    return events;
}

// TODO: This must not be done like this. It would take a gazillion
// operations to compare 2 keys of the map. Discuss the meaning of
// that later and come up with a better solution.
map<Events, Event> mapping_data_with_runs_events(int run, vector<Event> signal) {
    map<Events, Event> mapped_signal;
    for(size_t i = 0; i < signal.size(); i++) {
        Events run_event;
        run_event.run = run;
        run_event.event = i;
        mapped_signal.insert({run_event, signal[i]});
    }
    return mapped_signal;
}

int SignalProcessing_v2() {
    Channel chan = {
        .t = CHANNEL_SLOW,
        .n = 2,
    };
    vector<data_t> data = read_data_by_run(chan, RUN_START);
    vector<Event> events = split_data_to_events(chan, data);

    // TODO: Discuss this thingy
    auto ssignal = mapping_data_with_runs_events(RUN_START, events);
    for (const auto& pair : ssignal) {
        cout << pair.first.run << " " << pair.first.event << " " << pair.second[9999].y << endl;
    }
    return 0;
}

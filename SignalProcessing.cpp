#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <TLatex.h>

//----------------------//
using namespace std;
// ---------------------//

int run_start = 1;
int run_stop = run_start + 200;
int event_start = 0;
int event_stop = event_start + 14;
int threshold_slow = /*200*/ 40;
int threshold_fast = 20;

int one_alfa_l = 500E3;
int one_alfa_r = 800E3;
int one_peak_l = -10E3;
int one_peak_r = 10000E3;
int colib_l = 750E3;
int colib_r = 100000E3;
int left_line_integrate = 32000;
int right_line_integrate = 32700;
int seperate_coeff = 10;
bool colibration_mode = 0;
bool alfa = 0;					
string file_name = "D:\\Data\\old_setup\\2023\\230111_caen_archive\\f1";
string data_path = "C:\\Users\Mikheev\\Desktop\\code_root\\241112\\out_runNumb_eventNumb.txt";
FILE* outFile = NULL;
//----------------------//

const int points_per_event = 9999;
const int events_per_file = 1000;
const int sec_per_point = 16;
const int time_to_integrate = 32050;

//---------------------//

struct Point
{
	int x;
	int y;
};

struct Events 
{
	int runs;
	int events;
	
};

struct RunMode
{
	int start;
	int	stop;
};

struct data_bank_HIST 
{
	vector <double> base_line_avr;
	vector <double> base_line_sigma;
	vector <double> area;
	vector <double> x_max;
	vector <double> y_max; 	
	vector <double>	total_area;
	vector <Events> runs_events;
	vector <Events> runs_events_multy;
	
};

struct data_bank_SHORT_LIVE
{
	vector <double> area_colibr;
	vector <double> total_area_colibr;
	vector <double> x_max_colibr;
	vector <double> y_max_colibr;
	vector <double> area_short_live;
	vector <double> x_max_short_live;
	vector <double> y_max_short_live;
	vector <double> integrated_area;
	vector <double> integrated_area_piedistal;
	vector <double> integrated_area_no_piedistal;
	
};

struct data_bank_EVENT
{
	vector <vector<Point>> peaks;
	vector <vector<Point>> peaks_max;
	
};

data_bank_HIST HIST;
data_bank_EVENT PEAK;
data_bank_SHORT_LIVE COLIB;

void get_runs_events(int run, int event)
{
	Events s;
	s.runs = run;
	s.events = event;
	HIST.runs_events.push_back(s);
	
}

void get_runs_events_for_multy_signal(int _run, int _event)
{
	Events ss;
	ss.runs = _run;
	ss.events = _event;
	HIST.runs_events_multy.push_back(ss);
}

void get_base_line_data(double base_line_avr,double base_line_sigma)
{
	HIST.base_line_avr.push_back(base_line_avr);
	HIST.base_line_sigma.push_back(base_line_sigma);
	
}
void get_peaks(vector <vector<Point>> p,vector <vector<Point>> p_max)
{
	PEAK.peaks = p;
	PEAK.peaks_max = p_max;
	
}

void get_area_max(double area, double max_x, double max_y)
{
	HIST.area.push_back(area);
	COLIB.area_short_live.push_back(area);
	
	HIST.x_max.push_back(max_x);
	COLIB.x_max_short_live.push_back(max_x);
	
	HIST.y_max.push_back(max_y);
	COLIB.y_max_short_live.push_back(max_y);
	
}

void get_integrated_area(double area_integrated)
{
	COLIB.integrated_area.push_back(area_integrated);
	
}

void copy_data_colibr(void)
{
	 COLIB.area_colibr.insert(COLIB.area_colibr.end(), COLIB.area_short_live.begin(), COLIB.area_short_live.end());
	 COLIB.x_max_colibr.insert(COLIB.x_max_colibr.end(), COLIB.x_max_short_live.begin(), COLIB.x_max_short_live.end());
	 COLIB.y_max_colibr.insert(COLIB.y_max_colibr.end(), COLIB.y_max_short_live.begin(), COLIB.y_max_short_live.end());
} 

void get_tot_area_for_colibr(double total_area)
{
	COLIB.total_area_colibr.push_back(total_area);
}


void get_tot_area(double total_area)
{
	HIST.total_area.push_back(total_area);
}


vector <Events> read_data_file(string _data_path)
{
ifstream inputFile(_data_path);
vector <Events> v_adress;
string line;
	while(getline(inputFile, line))
	{
		istringstream iss(line);
		Events adress;
		if(iss >> adress.runs >> adress.events)
		{
			v_adress.push_back(adress);
		}
	}
	inputFile.close();
	
return v_adress;
} 

int columns (int number_events)
{
	int col;
	if (number_events == 1) col = 1;
	else col = 3;
	if (number_events > 9) col = 4;
	if (number_events > 12) col = 5;
	return col;
}

double* extract_xs(vector<Point> event)
{
	double* xs = new double[event.size()];
	for (int i = 0; i < event.size(); i++) 
	{
		xs[i] = sec_per_point*event[i].x;
	}
	return xs;
}

void clear_vector()
	{
		HIST.total_area.clear(); 
		HIST.total_area.shrink_to_fit();
		HIST.x_max.clear();
		HIST.x_max.shrink_to_fit();
		HIST.base_line_sigma.clear();
		HIST.base_line_sigma.shrink_to_fit();
		HIST.base_line_avr.clear();
		HIST.base_line_avr.shrink_to_fit();
		COLIB.total_area_colibr.clear();
		COLIB.total_area_colibr.shrink_to_fit();		
	}

void clear_vector_colib()
{
		COLIB.area_short_live.clear();
		COLIB.area_short_live.shrink_to_fit();
		COLIB.x_max_short_live.clear();
		COLIB.x_max_short_live.shrink_to_fit();
		COLIB.y_max_short_live.clear();
		COLIB.y_max_short_live.shrink_to_fit();	
}	

void clear_vector_area()
{
	HIST.area.clear();
	HIST.area.shrink_to_fit();
	COLIB.area_colibr.clear();
	COLIB.area_colibr.shrink_to_fit();
}

double* extract_ys(vector<Point> event)
{
	double* ys = new double[event.size()];
	for (int i = 0; i < event.size(); i++)
	{
		ys[i] = event[i].y;
	}
	return ys;
}

void alfa_mode(string mode)
{
	if (mode == "FAST")
	{
		alfa = 0;
		cout << "mode_cut is = " << alfa << endl;
	}
	else
	{
		cout << "alfa mode is = ? " << endl 
		<< " '1' - on " << endl << " '0' - off " << endl;
		cin >> alfa;
		cout << " mode_cut is = " << alfa << endl;	
	}
	
}
void colibration(string mode)
{
	if(mode == "FAST")
	{
		colibration_mode = 0;
		cout << "colibration_mode = " << colibration_mode << endl;
	}
	else 
	{
		cout << "colibration_mode = ?" << endl;
		cin >> colibration_mode;
		cout << "colibration_mode is = " << colibration_mode << endl;	
	}
	
}


RunMode conditions_mode(string choose_mode)
{
	RunMode Run;

	if (choose_mode == "HIST") 
	{
		Run.start = run_start;
		Run.stop = run_stop;	
	}
	if (choose_mode == "EVENT")
	{	
		Run.start = run_start;
		Run.stop  = run_start;
	}
	return Run;
}

//---------------------------------------//


vector <Point> read_data(int ch, int run, string choose_mode)
{	
	vector <Point> signal;
	signal.reserve((event_stop-event_start+1)*points_per_event);
	
	int coeff = 1;
	if(ch == 5||ch == 6) coeff = -1; 
	if(ch == 14||ch == 15) coeff = -1;
	if (ch >= 32)  coeff = -1;
	
	FILE *f = NULL;
	string filepath;
	
	ostringstream _ch;
	_ch << ch;
	string channel =  "__ch_" + _ch.str() + ".dat";

	ostringstream _run;
	_run << run; 
	filepath = file_name + "\\run_" + _run.str() + channel; 
	
	f = fopen(filepath.c_str(), "rb");
	
	if (!f) {
		cout << "Error opening file: " << filepath << endl;
	}
	
	int start; int stop;
	if (choose_mode == "EVENT")
	{
		fseek(f, 2 * event_start * points_per_event, SEEK_SET);
		start = event_start;
		stop = event_stop;		
	}
	if (choose_mode == "HIST")
	{
		start = 0;
		stop = events_per_file-1;
	}
	
	for (int event = start; event <= stop; event++) 
	{	
		vector<int16_t> buffer(points_per_event);
		size_t read_count = fread(buffer.data(), 2, points_per_event, f);			
		
		for (int j = 0; j < read_count; j++) 
		{ 
			Point point;	 	
			point.x = j;
			point.y = coeff*buffer[j];
			signal.push_back(point);
		}
	}
	
	fclose(f);
	cout << filepath << endl;
		
	return signal;
	
}

vector<Point> read_data_fast_tau(int run_vector, int ch, vector <Events> eventss)
{	
	vector <Point> signal;
	
	int coeff = 1;
	if(ch == 5||ch == 6) coeff = -1; 
	if(ch == 14||ch == 15) coeff = -1;
	if (ch >= 32)  coeff = -1;
	
	FILE *f = NULL;
	string filepath;
	
	ostringstream _ch;
	_ch << ch;
	string channel =  "__ch_" + _ch.str() + ".dat";
	
	//cout << "i = " << run_vector << endl; 
	ostringstream _run;
	_run << eventss[run_vector].runs; 
	filepath = file_name + "\\run_" + _run.str() + channel; 
	
	f = fopen(filepath.c_str(), "rb");
	
	if (!f) 
	{
		cout << "Error opening file: " << filepath << endl;
	}

	vector<int16_t> buffer(points_per_event);
	fseek(f, 2 * eventss[run_vector].events * points_per_event, SEEK_SET);
	size_t read_count = fread(buffer.data(), 2, points_per_event, f);			

	for (int j = 0; j < read_count; j++) 
	{ 
		Point point;	 	
		point.x = j;
		point.y = coeff*buffer[j];
		signal.push_back(point);
	}
	
	fclose(f);
	
	return signal;
}

vector<vector<Point>> partition_data(vector <Point> Signal, string choose_mode)
{
	vector <vector<Point>> signal_events;
	
	int size_vector;
	
	if (choose_mode == "EVENT")
	{
		size_vector = event_stop - event_start + 1;

		for (int event = 0; event < size_vector; event++)
		{
			vector <Point> temp(Signal.begin() + event*points_per_event, Signal.begin() + (event+1)*points_per_event);
			signal_events.push_back(std::move(temp));
		}	
	}
	if (choose_mode == "HIST")
	{
		size_vector = events_per_file;
		
		for (int event = 0; event < size_vector; event++)
		{
			vector <Point> temp(Signal.begin() + event*points_per_event, Signal.begin() + (event+1)*points_per_event);
			signal_events.push_back(std::move(temp));
		}
	}
	
	if (choose_mode == "SELECTION")
	{
		signal_events.push_back(Signal);
	}
	
	return signal_events;
}


vector <vector<Point>> normalize_baseLine(vector <vector<Point>> Signal, string choose_mode)
{	
	int period = 10000/sec_per_point;
	
	vector <vector<Point>> normalize_signal = Signal;
	double base_line_avr = 0;
	double base_line_sigma = 0;
	int sum = 0;
	int sum_sigma = 0;
	
	for (int event = 0; event < Signal.size(); event++)
	{
		sum = 0;
		sum_sigma = 0;
		
		for (int points = 0; points < period; points++)
		{
			sum += normalize_signal[event][points].y;
		}

		base_line_avr = sum/period;		
		
		for (int points = 0; points < period; points++) 
		{
			sum_sigma += pow(normalize_signal[event][points].y - base_line_avr, 2);
		}
		
		base_line_sigma = sqrt(sum_sigma/period);
		
		for (int points = 0; points < points_per_event; points++)
		{
			normalize_signal[event][points].y -= base_line_avr;
		}
		
		if (choose_mode == "HIST")
		{
			get_base_line_data(base_line_avr, base_line_sigma);
			
		}
	}	
	return normalize_signal;
}

//TODO вырезать альфа частицу только из ворот

double integrate_area_in_diapozone(int coeff_separation, vector<vector<Point>>& signal, int threshold, const int& event)
{
	int coeff = 1;
	double integrated_area = 0;
	
	for (int i = left_line_integrate/sec_per_point; i < right_line_integrate/sec_per_point; i++)
	{
		if (signal[event][i-1].y < threshold && signal[event][i].y > threshold)
		{
			coeff = coeff_separation;
			break;	
		}
		else 
			coeff = 1;
	}
	
	for (int j = left_line_integrate/sec_per_point; j <= right_line_integrate/sec_per_point; j++)
	{
		integrated_area += sec_per_point*(signal[event][j].y + signal[event][j-1].y)/2;
	}
	
	return coeff*integrated_area;
}

// void peak_fineder(vector <vector<Point>>& norm_signal, string choose_mode, int threshold, vector<Point>& peak, vector<Point>& peak_max, double& area, double& total_area, Point& max ,const int& event, const int& p)
// {
	// double prevmin = 0;
	// double prevmax = 0; 
			
	// if (norm_signal[event][p-1].y < threshold && norm_signal[event][p].y > threshold)
	// {
		// int p_min = p;
		
		// while ( norm_signal[event][p_min].y > threshold/6 && p_min > 0) p_min--; 

		// int p_max = p;
	
		// while ( norm_signal[event][p_max].y > threshold/5 && p_max < points_per_event) p_max++; 
		
		// if (p_min != prevmin && p_max != prevmax)
		// {
			// for (int i = p_min; i <= p_max; i++)
			// {
				// area += sec_per_point*(norm_signal[event][i].y + norm_signal[event][i-1].y)/2;
					
				// if (max.y < norm_signal[event][i].y)
					// max = norm_signal[event][i];
				
				// if (choose_mode == "EVENT")	
					// peak.push_back(norm_signal[event][i]);
				
			// }
		// }
		// else 
			// area = 0;
		
		// if (choose_mode == "EVENT")	
			// peak_max.push_back(max);
		
		// if (choose_mode == "HIST" && area > 0)
			// total_area += area;
								
		// if (choose_mode == "HIST" && area > 0)
		// {			
			// if ((max.y > threshold_slow + 20) && ((area > one_peak_l && area < one_peak_r) || !alfa))
				// get_area_max(area, sec_per_point * max.x, max.y);	
		// }	
		
		// prevmin = p_min;
		// prevmax = p_max;	
	// }

// }

void Find_peaks(vector <vector<Point>> norm_signal, string choose_mode, int run, string come_again, string integrate_status, int threshold)
{					
	double area = 0;
	double total_area = 0;
	double integrated_area = 0;
	double prevmin = 0;
	double prevmax = 0; 
	
	vector <Point> peak;
	vector <Point> peak_max;
	vector <vector<Point>> peakls;
	vector <vector<Point>> peakls_max;
	
	if (come_again == "SLOW")
	 outFile = fopen(data_path.c_str(), "a+");
	
	for (int event = 0; event < norm_signal.size(); event++)
	{
		
		if (come_again == "SLOW")
			integrated_area = integrate_area_in_diapozone(seperate_coeff, norm_signal, threshold, event);
	
		for (int p = 1; p < norm_signal[0].size(); p++)
		{
			Point max = {0,0};
			
		//	peak_fineder(norm_signal, choose_mode, threshold, peak, peak_max, area, total_area, max, event, p);	
			
			if (norm_signal[event][p-1].y < threshold && norm_signal[event][p].y > threshold)
			{
				int p_min = p;
				
				while ( norm_signal[event][p_min].y > threshold/6 && p_min > 0) p_min--; 

				int p_max = p;
			
				while ( norm_signal[event][p_max].y > threshold/5 && p_max < points_per_event) p_max++; 
				
				if (p_min != prevmin && p_max != prevmax)
				{
					for (int i = p_min; i <= p_max; i++)
					{
						area += sec_per_point*(norm_signal[event][i].y + norm_signal[event][i-1].y)/2;
							
						if (max.y < norm_signal[event][i].y)
							max = norm_signal[event][i];
						
						if (choose_mode == "EVENT")	
							peak.push_back(norm_signal[event][i]);
						
					}
				}
				else 
					area = 0;
				
				if (choose_mode == "EVENT")	
					peak_max.push_back(max);
				
				if (choose_mode == "HIST" && area > 0)
					total_area += area;
										
				if (choose_mode == "HIST" && area > 0)
				{			
					if ((max.y > threshold_slow + 20) && ((area > one_peak_l && area < one_peak_r) || !alfa))
						get_area_max(area, sec_per_point * max.x, max.y);	
				}	
				
				prevmin = p_min;
				prevmax = p_max;	
			}
			
			if (choose_mode == "HIST")
			{
				if ((max.x > (left_line_integrate)/sec_per_point && max.x < (right_line_integrate)/sec_per_point))
				{
					if (!colibration_mode || (total_area > 0 && (total_area > colib_r || total_area < colib_l)))
					{
						get_runs_events_for_multy_signal(run, event);							
					}
						
				}
				
			}
			
			
		
		area = 0;
			
		} //цикл точек 
		
		bool event_shoud_be_cut;
		
		if (come_again == "SLOW")
		{
			
			event_shoud_be_cut = 1;
			for (int i = 0; i < COLIB.area_short_live.size(); i++)
			{
				if (COLIB.area_short_live[i] > one_peak_r && (COLIB.x_max_short_live[i] > left_line_integrate && COLIB.x_max_short_live[i] < right_line_integrate))
				{
					event_shoud_be_cut = 0;
				//	cout << event << "\t" << " area = " << COLIB.area_short_live[i] << "\t" << " x = " << COLIB.x_max_short_live[i] << endl;
					break;
				}
				
				//cout << " max.x_origin_vector = " << COLIB.area_short_live[i] << endl;	
				
			}
			
			//if(event_shoud_be_cut || !colibration_mode)
			
		}
		
		
		if (choose_mode == "EVENT")
		{	 
			peakls.push_back(peak);
			peakls_max.push_back(peak_max);
			
			get_peaks(peakls, peakls_max);
			
			peak.clear();
			peak_max.clear();
			peak.shrink_to_fit();
			peak_max.shrink_to_fit();
			
		}
	
		if (choose_mode == "HIST" && total_area > 0)//&& total_area > 0
		{	
			if (total_area > one_alfa_l && total_area < one_alfa_r || !alfa)
			{
				get_tot_area(total_area);
			
				if (come_again == "SLOW")
				{
					ostringstream output_data;
					
					output_data << run << "\t" << event << "\n";
					fwrite(output_data.str().c_str(), 1, output_data.str().length(), outFile);
					
					get_runs_events(run, event);						
				}
			} 
		}
		
		if (!colibration_mode || (total_area > 0 && (total_area > colib_r || total_area < colib_l)))
		{
			get_tot_area_for_colibr(total_area);
			
			if (come_again == "SLOW")
			{
				copy_data_colibr();
				
				//if(event_shoud_be_cut)	
					get_integrated_area(integrated_area);			
				
			}
			
		}
		
		clear_vector_colib();
		total_area = 0;
		
		
	}//цикл событий
	
	if (come_again == "SLOW")
	{
		fclose(outFile);
	}
	
}

void Plot_EVENT(vector<vector<Point>> divided_data, int columns, int threshold)
{
	TCanvas *c2 = new TCanvas("c2","multipads",1000,1000);
	
	int rows = ceil((float)divided_data.size() / columns);
	c2->Divide(columns,rows,0.01,0.01);
	
	
	for(int i = 0; i < divided_data.size(); i++)
	{
		c2->cd(i+1);
		vector<Point> event = divided_data[i];
		double* x_signal = extract_xs(event);
		double* y_signal = extract_ys(event);
		if (event.size() > 0) 
		{
		TGraph *signal = new TGraph(event.size(), x_signal, y_signal);
		signal->SetMarkerStyle(10);
		signal->Draw();
		ostringstream title;
		title << "Event: " << event_start + i; 
		signal->SetTitle(title.str().c_str());
		signal->GetXaxis()->SetTitle("Time");
		signal->GetYaxis()->SetTitle("signal level");
		
		
		TF1 *line = new TF1("line","[0] + [1]*x", x_signal[0], x_signal[event.size()-1]);
		line->SetParameter(0, threshold);
		line->SetParameter(1, 0);
		line->SetNpx(10000);
		line->Draw("L same");
		}
		
		vector <Point> event_peaks = PEAK.peaks[i];
		if (PEAK.peaks.size() > 0) 
		{
			TGraph*signal_peaks = new TGraph(PEAK.peaks[i].size(), extract_xs(event_peaks),  extract_ys(event_peaks));
			signal_peaks->Draw("P same");
			signal_peaks->SetMarkerColor(kOrange+1);
			signal_peaks->SetMarkerStyle(20);
		}
		
		vector <Point> event_peaks_max = PEAK.peaks_max[i]; 
		if (event_peaks_max.size() > 0) 
		{
			TGraph*signal_peaks_max = new TGraph(PEAK.peaks_max[i].size(), 
			extract_xs(event_peaks_max),  extract_ys(event_peaks_max));
			
			signal_peaks_max->Draw("P same");
			signal_peaks_max->SetMarkerColor(kBlue+2);
			signal_peaks_max->SetMarkerStyle(20);
		}
		
	}
	
}

void Plot_HIST(const char* canvas_name)
{
	int rangeSignal = 4E6;
	if (!(strcasecmp(canvas_name, "FAST_signal"))) 	
	{
		rangeSignal = 180E3;
	}
	TH1* h1_base_line_sigma = new TH1F("h1_base_line_sigma", "h1_base_line_sigma", 1000, -50, 50);	     
	TH1* h1_base_line_avr = new TH1F("h1_base_line_avr", "h1_base_line_avr", 1000, -5000, 5000);  
	TH1F* h1_peak_time_weight = new TH1F("h1","h1_avr_weight",1600,0,160000);
	TH1F* h1_peak_area = new TH1F("h1", "h1_peak_area", 550*5, -100000, 1.5E6);  
	TH1F* h1_total_area = new TH1F("h1", "h1_total_area", 2000*5, -200000, 10E6); 
	TH2F* hist_peak_area_ev_vs_evnum = new TH2F("h2", "h2", 110, 0, HIST.total_area.size(), 100, 0, rangeSignal);
	
	TCanvas *c1 = new TCanvas(canvas_name,canvas_name,1000,1000);
	c1->Divide(3,2,0.01,0.01);
	
	string str(canvas_name);
	
	
	
	for (int run = 0; run < HIST.total_area.size(); run++)
	{
		hist_peak_area_ev_vs_evnum->Fill(run, HIST.total_area[run]);	
	}
	
	c1->cd(1);
	c1->cd(1)->SetTickx();
	c1->cd(1)->SetTicky();
	hist_peak_area_ev_vs_evnum->GetXaxis()->SetTitle("event");
	hist_peak_area_ev_vs_evnum->GetYaxis()->SetTitle("peak_area_ev [ADC_code*ns]");
	hist_peak_area_ev_vs_evnum->Draw("colz");
	
	TProfile* prof_hist_peak_area_ev_vs_evnum = hist_peak_area_ev_vs_evnum->ProfileX();	
	//prof_hist_peak_area_ev_vs_evnum->Draw("same");
	prof_hist_peak_area_ev_vs_evnum->SetMarkerStyle(20);
	prof_hist_peak_area_ev_vs_evnum->SetMarkerColor(kRed);		
	


	c1->cd(2);
	c1->cd(2)->SetTickx();
	c1->cd(2)->SetTicky();
	
	cout << "size_total_colibr = " << COLIB.total_area_colibr.size() << endl;
	cout << "size_area_colibr = " << COLIB.area_colibr.size() << endl;
	cout << "size_area = " << HIST.area.size() << endl;
	cout << "size_xmax = " << HIST.x_max.size() << endl;
	
	if (!colibration_mode) 
	{
		for (int i = 0; i < HIST.total_area.size(); i++)
		{
			h1_total_area->Fill(HIST.total_area[i]);	
		}
			
	}
	else
	{
		for (int i = 0; i < COLIB.total_area_colibr.size(); i++)
		{
			h1_total_area->Fill(COLIB.total_area_colibr[i]);	
		}
		
	} 
	
	h1_total_area->GetXaxis()->SetTitle("Area");
	h1_total_area->GetYaxis()->SetTitle("Number");
	h1_total_area->SetFillColor(kBlue-1.5);
	h1_total_area->Draw();
	// TF1*fit_gaus_0 = new TF1("fit_gaus_0", "gaus", one_alfa_l, one_alfa_r);
	// h1_total_area->Fit(fit_gaus_0, "R");
	// double area_alfa_area = fit_gaus_0->Integral(one_alfa_l, one_alfa_r);
	// cout << " area_alfa_area = " << area_alfa_area << endl;
	
	c1->cd(3);
	c1->cd(3)->SetTickx();
	c1->cd(3)->SetTicky();
	
	if (!colibration_mode)
	{
		for (int i = 0; i < HIST.area.size(); i++)
		{
			h1_peak_area->Fill(HIST.area[i]);	
		}
	}
	else 
	{
		for (int i = 0; i < COLIB.area_colibr.size(); i++)
		{
			h1_peak_area->Fill(COLIB.area_colibr[i]);	
		}
	}
	
	h1_peak_area->GetXaxis()->SetTitle("Area");
	h1_peak_area->GetYaxis()->SetTitle("Number of peaks");
	h1_peak_area->SetFillColor(kBlue-1.5);
	h1_peak_area->Draw();
	// TF1*fit_landau = new TF1("fit_landau", "landau", one_peak_l, one_alfa_r);
	// h1_peak_area->Fit(fit_landau);
	// double area_peak_area = fit_landau->
	// Integral(one_peak_l, one_alfa_r);
	// double mean_area_peak_area = area_peak_area / abs(one_alfa_r - one_peak_l);
	// cout << " mean_area_peak_area = " << mean_area_peak_area << "\t" << "mean = " << fit_landau->GetParameter(1) << endl;
	
	if (!(strcasecmp(canvas_name, "SLOW_signal")))
	{
		ofstream outFile("C:\\Users\\Mikheev\\Desktop\\code_root\\phe_Hist.txt");
		if (outFile.is_open())
		{
			cout << "file phe_Hist.txt is open" << endl;
			
			for (int bin = 1; bin <= h1_peak_area->GetNbinsX(); bin++) 
			{ 
				double binContent = h1_peak_area->GetBinContent(bin);
				int intBinContent = static_cast<int>(binContent); 
				double binCenter = h1_peak_area->GetBinCenter(bin);
				outFile << static_cast<int>(binCenter) << "\t" << intBinContent << endl;
			
			}
		}
		outFile.close();	
		cout << "data save complete " << endl;
		
	}
	
		
	c1->cd(4);
	c1->cd(4)->SetTickx();
	c1->cd(4)->SetTicky();
	
	if (!colibration_mode)
	{
		for (int i = 0; i < HIST.area.size(); i++) 
		{
			h1_peak_time_weight->Fill(HIST.x_max[i], HIST.area[i]);
		}
	}
	else 
	{
		for (int i = 0; i < COLIB.area_colibr.size(); i++) // как сделать ?
		{
			h1_peak_time_weight->Fill(COLIB.x_max_colibr[i], COLIB.area_colibr[i]);
		}
		
	}
	
	
	h1_peak_time_weight->Draw("E");
	gPad->SetLogy();
	
	c1->cd(5);
	c1->cd(5)->SetTickx();
	c1->cd(5)->SetTicky();
	for (int i = 0; i < HIST.base_line_sigma.size(); i++)
	{
		h1_base_line_sigma->Fill(HIST.base_line_sigma[i]);
	}
	h1_base_line_sigma->SetFillColor(kBlue-1.5);
	h1_base_line_sigma->Draw();
	gPad->SetLogy();
	
	c1->cd(6);
	c1->cd(6)->SetTickx();
	c1->cd(6)->SetTicky();
	for (int i = 0; i < HIST.base_line_avr.size(); i++)
	{
		h1_base_line_avr->Fill(HIST.base_line_avr[i]);
	}
	h1_base_line_avr->Draw();
	gPad->SetLogy();

}

void Plot_multy_signals(vector <vector<Point>> normalize_signal_multy, TCanvas* canvas, int left_line_integrate, int left_right_integrate)
{
	canvas->cd(3);
	canvas->cd(3)->SetTickx();
	canvas->cd(3)->SetTicky();
	
	int roof_signal = 1800;
	vector <Point> events;
	
	int start_events_runs = 40;
	int number_of_signals = 20;
	
	int base_kolor = 0;
	int counter_for_color = 0;
	int mix = 0;
	double color_plus = 0;
	
	for (int i = start_events_runs; i <= start_events_runs + number_of_signals; i++) /*_normalize_signal_multy.size()/points_per_event*/
	{
		events = normalize_signal_multy[i];
		TGraph *signal = new TGraph(events.size(), extract_xs(events), extract_ys(events));
		signal->SetMarkerStyle(15);
	
		// counter_for_color++;
		// if (i == 1) base_kolor = 632; 
		// if (i == 6) base_kolor = 416; 
		// if (i == 4) base_kolor = 600;
		// if (i == 3)	base_kolor = 800;
		// if (i == 5)	base_kolor = 616;
		// if (i == 2)	base_kolor = 900;
		// signal->SetLineColor(base_kolor);
		
		int color = TColor::GetColor(50 + (i * 50) % 200, 50 + (i * 100) % 200, 50 + (i * 150) % 200); // Генерация цвета
		signal->SetLineColor(color);
		signal->SetMarkerColor(color);
		
		if (i == start_events_runs) signal->Draw();
		if (i > start_events_runs) signal->Draw("SAME");
		
		signal->GetXaxis()->SetTitle("Time");
		signal->GetYaxis()->SetTitle("signal level");
		signal->GetYaxis()->SetRangeUser(0, roof_signal); // Устанавливаем границы по оси Y	
		
		//if (counter_for_color == 3) counter_for_color = 0; mix += 1;
	}
	
	TLine *line_lefty = new TLine(left_line_integrate, 0, left_line_integrate, roof_signal);
	line_lefty->SetLineColor(kRed);
	line_lefty->Draw("same");
	line_lefty->SetLineStyle(2);
	line_lefty->SetLineWidth(2);
	
	TLine *line_rigty = new TLine(right_line_integrate, 0, right_line_integrate, roof_signal);
	line_rigty->SetLineColor(kRed);
	line_rigty->Draw("same");
	line_rigty->SetLineStyle(2);
	line_rigty->SetLineWidth(2);
	
	
}

// void Plot_intagrated_Area(vector <vector<Point>> colected_signal, TCanvas* canvas)
// {
	// for (int i = 0; i < colected_signal.size(); i++)
	// {
		// for ()
	// }
	
// }


void Multy_Signals(TCanvas* canvas, int ch)
{
	vector <Point> multy_signal;
	vector <vector<Point>> partition_signal;
	vector <vector<Point>> normalize_signal;
	vector <vector<Point>> colect_multy_signal;
	
	cout << "HIST.runs_events_multy = " << HIST.runs_events_multy.size() << endl;

	for (int i = 0; i < HIST.runs_events_multy.size(); i++)
	{
		multy_signal = read_data_fast_tau(i, ch, HIST.runs_events_multy);
		partition_signal = partition_data(multy_signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "NO_GIST");	
		colect_multy_signal.push_back(normalize_signal[0]);
	}
	Plot_multy_signals(colect_multy_signal, canvas, left_line_integrate, right_line_integrate);

}

void Plot_area_Amp(TCanvas* canvas)
{
	
	if (!colibration_mode)
	{
		
		TH2F* AMP_AREA = new TH2F("h2", "AMP_AREA", 100, 0, 150000, 100, 0, 700);	
		
		canvas->cd(1);
		canvas->cd(1)->SetTickx();
		canvas->cd(1)->SetTicky();
		
		for (int i = 0; i < HIST.area.size(); i++) 
		{
			AMP_AREA->Fill( HIST.area[i], HIST.y_max[i]);
		}
		AMP_AREA->GetXaxis()->SetTitle("area_of_phe");
		AMP_AREA->GetYaxis()->SetTitle("Amplitude_of_peak");
		AMP_AREA->Draw("COLZ");
	}
	else 
	{
		TH2F* AMP_AREA = new TH2F("h2", "AMP_AREA", 100, 0, 150000 , 100, 0, 700);
		
		canvas->cd(1);
		canvas->cd(1)->SetTickx();
		canvas->cd(1)->SetTicky();
		
		for (int i = 0; i < COLIB.area_colibr.size(); i++) 
		{
			AMP_AREA->Fill(COLIB.area_colibr[i], COLIB.y_max_colibr[i]);
		}
		
		AMP_AREA->GetXaxis()->SetTitle("area_of_phe");
		AMP_AREA->GetYaxis()->SetTitle("Amplitude_of_peak");
		AMP_AREA->Draw("COLZ");
	}
	

}


void Plot_intagrate_Area(TCanvas* canvas, int seperate_coeff)
{
	TH1F* h1_hist_inegrate_area = new TH1F("h3", "h1_hist_inegrate_area", 4088, -1000000, 7E6); 
	
	canvas->cd(2);
	for (int i = 0; i < COLIB.integrated_area.size(); i++)
	{
		h1_hist_inegrate_area->Fill(COLIB.integrated_area[i]);
	}
	
	int seperate_line = 22000;																		//
	int lowbound = 0;
	int upbound = h1_hist_inegrate_area->FindBin(seperate_line);
	
	for (int i = 1; i <= h1_hist_inegrate_area->GetNbinsX(); ++i) {
        if (h1_hist_inegrate_area->GetBinContent(i) > 0) {
            lowbound = i;
            break;
        }
    }
	
	double Integral_all = h1_hist_inegrate_area->Integral();
	double Integral_piedistal = h1_hist_inegrate_area->Integral(lowbound, upbound);
	double Probability_of_Piedistal = Integral_piedistal/Integral_all;
	double mean_S = 0;
	double S_piedistal = 0;
	double S_signal = 0;
	double total_S = 0;
	double mean_phe = 0;
	double G = 0;
	int total_Events = 0;
	int total_Events_signal = 0;
	
	for (int i = h1_hist_inegrate_area->FindBin(-20000); i < h1_hist_inegrate_area->FindBin(seperate_line); i++)
	{
		S_piedistal += h1_hist_inegrate_area->GetBinContent(i) * h1_hist_inegrate_area->GetBinCenter(i);
		total_Events += h1_hist_inegrate_area->GetBinContent(i);
	}
	
	for (int i = h1_hist_inegrate_area->FindBin(seperate_line); i < h1_hist_inegrate_area->GetNbinsX(); i++)
	{
		S_signal += h1_hist_inegrate_area->GetBinContent(i) * h1_hist_inegrate_area->GetBinCenter(i);
		total_Events_signal += h1_hist_inegrate_area->GetBinContent(i);
	}
	
	S_piedistal = S_piedistal/total_Events;
	S_signal = S_signal/total_Events_signal;
	S_signal = S_signal/seperate_coeff;
	
	total_S = (total_Events*S_piedistal + total_Events_signal*S_signal)/(total_Events + total_Events_signal);
	mean_S = h1_hist_inegrate_area-> GetMean();
	mean_phe = -log(Probability_of_Piedistal);
	G = total_S/mean_phe; 	
	
	cout << "Integral_all = " << Integral_all << endl;
	cout << "Integral_piedistal = " << Integral_piedistal << endl;
	cout << "Probability_of_Piedistal = " << Probability_of_Piedistal << endl;
	cout << "mean_S = " << mean_S <<endl;
	cout << "piedistal_S = " << S_piedistal <<endl;
	cout << "signal_S = " << S_signal <<endl;
	cout << "total_S = " << total_S <<endl;
	cout << "mean_phe = " << mean_phe <<endl;
	cout << "G = " << G  <<endl;


	TLatex latex;
	latex.SetTextSize(1); 
	latex.DrawLatex(1000, 200, to_string(static_cast<double>(G)).c_str());
	
	h1_hist_inegrate_area->GetXaxis()->SetTitle("area_of_phe");
	h1_hist_inegrate_area->GetYaxis()->SetTitle("event");
	h1_hist_inegrate_area->Draw();
	
	TLine *line_seperate = new TLine(seperate_line, 0 , seperate_line, h1_hist_inegrate_area->GetMaximum());
	line_seperate->Draw();
	line_seperate->SetLineColor(kRed);

	
}

void Other_PARAMETERS(int ch)
{
	TCanvas *c4 = new TCanvas("OTHER PARAMETERS","OTHER PARAMETERS",1000,1000);
	c4 ->Divide(3,2,0.01,0.01);
	
	Plot_area_Amp(c4);
	Plot_intagrate_Area(c4, seperate_coeff);
	Multy_Signals(c4, ch);
	
}

vector<Point> cut_normalize_signal(vector <vector<Point>> normalize_signal, int left_border, int right_border)
{
	vector<Point> cut_norm_signal;
	
	for (int i = 0; i < normalize_signal[0].size(); i++)
	{
		if (i >= left_border/sec_per_point && i <= right_border/sec_per_point)
		{
			cut_norm_signal.push_back(normalize_signal[0][i]);
		}
	}
	return cut_norm_signal;
}

void func_HIST(string HIST, int ch)
{
	colibration("SLOW");
	alfa_mode("SLOW");
	vector <Point> Signal;
	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;
	vector<vector<Point>> cutty_norm_signal;

	RunMode Run = conditions_mode(HIST);
	for (int run = Run.start; run <= Run.stop; run++)
	{ 
		Signal = read_data(ch, run, HIST);
		partition_signal = partition_data(Signal, HIST);
		normalize_signal = normalize_baseLine(partition_signal, HIST);
		cutty_norm_signal.push_back(cut_normalize_signal(normalize_signal, left_line_integrate, right_line_integrate));
		Find_peaks(normalize_signal, HIST, run, "SLOW", "dontINTEGRATE", threshold_slow);
		Find_peaks(cutty_norm_signal, HIST, run, "0","INTEGRATE", threshold_slow);
	}
	
	Plot_HIST("SLOW_signal");
	
	clear_vector();
}

void func_EVENT(string EVENT, int ch)
{
	int threshold_temp = threshold_slow;
	if (ch == 5 || ch == 6) threshold_temp = threshold_fast;
	if (ch == 2 || ch == 3) threshold_temp = threshold_slow; 
	
	colibration("FAST");
	alfa_mode("FAST");
	vector <Point> Signal;
	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;
	
	RunMode Run = conditions_mode(EVENT);
	for (int run = Run.start; run <= Run.stop; run++)
	{ 
		Signal = read_data(ch, run, EVENT);
		partition_signal = partition_data(Signal, EVENT);
		normalize_signal = normalize_baseLine(partition_signal, EVENT);
		Find_peaks(normalize_signal, EVENT, run, "0", "dontINTEGRATE", threshold_temp); 
		Plot_EVENT(normalize_signal, columns(event_stop-event_start+1), threshold_temp);
	}
}


void func_HIST_fast_tau(int ch)
{
	colibration("FAST");
	alfa_mode("FAST");
	vector<Point> signal; 
	vector <vector<Point>> partition_signal;
	vector <vector<Point>> normalize_signal;
	
	
	for (int i = 0; i < HIST.runs_events.size(); i++)
	{
		signal = read_data_fast_tau(i, ch, HIST.runs_events);
		partition_signal = partition_data(signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "HIST");
		Find_peaks(normalize_signal, "HIST", 0, "0", "dontINTEGRATE", threshold_fast);
	}
	Plot_HIST("FAST_signal");
	
	clear_vector();
}


int SignalProcessing() 
{
	outFile = fopen(data_path.c_str(), "w");
	fclose(outFile);
	
	clock_t start_first = clock();
	
	func_HIST("HIST", 1);
	
//	Other_PARAMETERS(1);
	
	clear_vector_area();
	
//	func_EVENT("EVENT", 14);
	
	func_HIST_fast_tau(14);	
	
	
	
	clock_t end_fisrt = clock();
	
	
	 
	 double duration = double(end_fisrt - start_first);
	 
	
	 cout << "time_all = "<< duration << endl;
	 
	 return 0;
}
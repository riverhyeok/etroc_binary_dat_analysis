#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <cmath>

// ROOT Headers for Histogram & PNG Export
#include <TROOT.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>

using namespace std;

// =====================================================================
// 1. ETROC2 Bit Masking & Physical Constants
// =====================================================================
const double T3_NS = 3.125; // TDC Reference Strobe 간격 (ns)

namespace etroc {
    static const uint64_t SYNC = 0x3C5C;
    inline bool is_data(uint64_t d)    { return  (d >> 39) & 1; }
    inline bool is_sync(uint64_t d)    { return !is_data(d) && (((d >> 24) & 0x7FFF) == SYNC); }
    inline bool is_header(uint64_t d)  { return  is_sync(d) && (((d >> 22) & 0x3) == 0x0); }
    inline bool is_filler(uint64_t d)  { return  is_sync(d) && (((d >> 22) & 0x3) == 0x2); }
    inline bool is_trailer(uint64_t d) { return !is_data(d) && !is_sync(d); }

    inline uint32_t hdr_bcid(uint64_t d)   { return  d        & 0xFFF; }
    inline uint8_t  hdr_type(uint64_t d)   { return (d >> 12) & 0x3;   }
    inline uint8_t  hdr_l1_cnt(uint64_t d) { return (d >> 14) & 0xFF;  }

    // Regular Data 비트 마스킹
    inline uint8_t  dat_ea(uint64_t d)     { return (d >> 37) & 0x3;   } 
    inline uint8_t  dat_col(uint64_t d)    { return (d >> 33) & 0xF;   }
    inline uint8_t  dat_row(uint64_t d)    { return (d >> 29) & 0xF;   }
    inline uint16_t dat_toa(uint64_t d)    { return (d >> 19) & 0x3FF; } // 10-bit TOA
    inline uint16_t dat_tot(uint64_t d)    { return (d >> 10) & 0x1FF; } // 9-bit TOT
    inline uint16_t dat_cal(uint64_t d)    { return d & 0x3FF;         } // 10-bit CAL

    inline uint8_t  trl_status(uint64_t d) { return (d >> 16) & 0x3F;   }
    inline uint16_t trl_hits(uint64_t d)   { return (d >>  8) & 0xFF;   }
    inline uint32_t trl_chipid(uint64_t d) { return (d >> 22) & 0x1FFFF;}

    uint8_t calculate_crc8(const std::vector<uint8_t>& data_bytes) {
        uint8_t crc = 0x00;
        for (uint8_t b : data_bytes) {
            crc ^= b;
            for (int i = 0; i < 8; i++) {
                if (crc & 0x80) crc = (crc << 1) ^ 0x2F;
                else            crc <<= 1;
            }
        }
        return crc;
    }
}

enum FrameType { HEADER = 0, DATA = 1, TRAILER = 2, FILLER = 3, UNKNOWN = 4 };

FrameType get_payload_type(uint64_t w) {
    if (etroc::is_data(w))    return DATA;
    if (etroc::is_header(w))  return HEADER;
    if (etroc::is_filler(w))  return FILLER;
    if (etroc::is_trailer(w)) return TRAILER;
    return UNKNOWN;
}

struct SystemStats {
    uint64_t total_32b_words = 0;
    uint64_t total_ram_buffer_empty = 0;
    void add(const SystemStats& o) {
        total_32b_words += o.total_32b_words;
        total_ram_buffer_empty += o.total_ram_buffer_empty;
    }
};

struct RecoStats {
    uint64_t total_frames = 0, headers = 0, trailers = 0, fillers = 0, hits_raw = 0, valid_hits = 0;
    uint64_t events_missing_header = 0, orphan_data_words = 0, words_dropped_at_start = 0;
    
    uint64_t wp_boundary_cuts = 0;
    uint64_t weird_wp_boundary_cuts = 0;
    
    uint64_t events_perfect_hits = 0, events_hit_mismatch = 0, events_corrupted_bcid = 0, events_discarded_by_jump = 0;
    uint64_t l1_counter_jumps = 0, l1_jump_sum = 0, match_intra_colrow = 0, drop_intra_colrow = 0;
    uint64_t match_hdr_bcid = 0, drop_hdr_bcid = 0, drop_hdr_bcid_only = 0, counter_a_jumps = 0, counter_a_jump_sum = 0;
    uint64_t ea_counts[4] = {0}, fatal_ea_errors = 0, crc_match = 0, crc_mismatch = 0;
    uint64_t expected_total_hits = 0, sw_lost_hits = 0;
    uint64_t header_type_counts[4] = {0}, buffer_status_counts[4] = {0};
    uint64_t hit_map[16][16] = {0};
    uint64_t fake_trailers_blocked = 0, unknown_frames = 0;
    uint64_t crc_err_ea = 0, crc_err_missing = 0, crc_err_missing_with_ea = 0, crc_err_silent = 0;
    uint64_t events_with_missing_hits = 0, events_missing_hits_and_crc_err = 0;
    uint64_t events_missing_hits_but_crc_ok = 0, events_with_ea_errors = 0, events_with_ea_and_crc_err = 0;

    void add(const RecoStats& o) {
        total_frames += o.total_frames; headers += o.headers; trailers += o.trailers; fillers += o.fillers;
        hits_raw += o.hits_raw; valid_hits += o.valid_hits; events_missing_header += o.events_missing_header;
        orphan_data_words += o.orphan_data_words; words_dropped_at_start += o.words_dropped_at_start;
        wp_boundary_cuts += o.wp_boundary_cuts; weird_wp_boundary_cuts += o.weird_wp_boundary_cuts;
        events_perfect_hits += o.events_perfect_hits; events_hit_mismatch += o.events_hit_mismatch;
        events_corrupted_bcid += o.events_corrupted_bcid; events_discarded_by_jump += o.events_discarded_by_jump;
        l1_counter_jumps += o.l1_counter_jumps; l1_jump_sum += o.l1_jump_sum;
        match_intra_colrow += o.match_intra_colrow; drop_intra_colrow += o.drop_intra_colrow;
        match_hdr_bcid += o.match_hdr_bcid; drop_hdr_bcid += o.drop_hdr_bcid; drop_hdr_bcid_only += o.drop_hdr_bcid_only;
        counter_a_jumps += o.counter_a_jumps; counter_a_jump_sum += o.counter_a_jump_sum;
        for(int i=0; i<4; i++) { ea_counts[i] += o.ea_counts[i]; header_type_counts[i] += o.header_type_counts[i]; buffer_status_counts[i] += o.buffer_status_counts[i]; }
        fatal_ea_errors += o.fatal_ea_errors; crc_match += o.crc_match; crc_mismatch += o.crc_mismatch;
        expected_total_hits += o.expected_total_hits; sw_lost_hits += o.sw_lost_hits;
        for(int r=0; r<16; r++) for(int c=0; c<16; c++) hit_map[r][c] += o.hit_map[r][c];
        fake_trailers_blocked += o.fake_trailers_blocked; unknown_frames += o.unknown_frames;
        crc_err_ea += o.crc_err_ea; crc_err_missing += o.crc_err_missing; crc_err_missing_with_ea += o.crc_err_missing_with_ea;
        crc_err_silent += o.crc_err_silent; events_with_missing_hits += o.events_with_missing_hits;
        events_missing_hits_and_crc_err += o.events_missing_hits_and_crc_err; events_missing_hits_but_crc_ok += o.events_missing_hits_but_crc_ok;
        events_with_ea_errors += o.events_with_ea_errors; events_with_ea_and_crc_err += o.events_with_ea_and_crc_err;
    }
};

// Global 데이터 보관소
std::map<uint8_t, TH1D*> g_toa; 
std::map<uint8_t, TH1D*> g_tot; 
std::map<uint8_t, TH1D*> g_cal; 
std::map<uint8_t, RecoStats> g_stats;
std::map<uint8_t, SystemStats> g_sys;
std::map<uint8_t, uint64_t> g_valid_hits;

// =====================================================================
// 2. Analyzer Class
// =====================================================================
class EtrocChannelAnalyzer {
    string name_, file_prefix_;
    uint32_t expected_chip_id_;
    RecoStats s_;
    uint8_t chan_tag_;
    
    TH1D* h_toa_;
    TH1D* h_tot_;
    TH1D* h_cal_;

    bool is_synced_ = false, expect_l1_jump_ = false, in_event_ = false, has_pending_event_ = false;
    uint64_t pending_trailer_w_ = 0, event_id_ = 0;
    uint32_t hdr_bcid_ = 0;
    uint8_t  hdr_type_ = 0, last_l1_cnt_ = 0;
    bool     is_first_event_ = true;
    
    std::vector<uint64_t> current_event_data_;
    std::vector<uint8_t> event_bytes_;
    std::vector<string> hit_stream_buffer_; 
    
    std::vector<uint64_t> recent_words_;
    int cut_print_count_ = 0;
    uint64_t global_word_index_ = 0; 

    void push_word_bytes(std::vector<uint8_t>& target_vec, uint64_t w) {
        for (int i = 4; i >= 0; i--) target_vec.push_back((w >> (i * 8)) & 0xFF);
    }

    void commit_pending_event() {
        s_.header_type_counts[hdr_type_]++;
        uint16_t trl_h = etroc::trl_hits(pending_trailer_w_);
        uint8_t status = etroc::trl_status(pending_trailer_w_);
        s_.buffer_status_counts[(status >> 4) & 0x3]++;
        uint64_t event_hits = current_event_data_.size();
        if (trl_h == 0 && event_hits >= 200) trl_h = 256;

        bool current_event_has_ea_error = false;

        for (uint64_t hit_w : current_event_data_) {
            uint8_t ea = etroc::dat_ea(hit_w);
            s_.ea_counts[ea & 0x3]++;
            
            if (ea >= 0x2) { 
                s_.fatal_ea_errors++; 
                current_event_has_ea_error = true; 
                continue; 
            }

            uint8_t col = etroc::dat_col(hit_w);
            uint8_t row = etroc::dat_row(hit_w);
            
            s_.valid_hits++;
            s_.hit_map[row][col]++;

            uint16_t tot_code = etroc::dat_tot(hit_w);
            uint16_t toa_code = etroc::dat_toa(hit_w);
            uint16_t cal_code = etroc::dat_cal(hit_w);

            double toa_ns = 0.0;
            double tot_ns = 0.0;
            
            if (cal_code > 0) {
                double t_bin = T3_NS / (double)cal_code;
                toa_ns = t_bin * (double)toa_code;
                
                double tot_corrected = (2.0 * tot_code) - std::floor(tot_code / 32.0);
                tot_ns = t_bin * tot_corrected;

                // 물리 단위 플롯 내부 메모리 적재
                h_toa_->Fill(toa_ns);
                h_tot_->Fill(tot_ns);
                g_toa[chan_tag_]->Fill(toa_ns);
                g_tot[chan_tag_]->Fill(tot_ns);
            }
            
            h_cal_->Fill(cal_code);
            g_cal[chan_tag_]->Fill(cal_code);
            g_valid_hits[chan_tag_]++;
            
            // 파이썬 시각화용 데이터에는 변환된 물리 시간 전달
            hit_stream_buffer_.push_back(
                to_string(event_id_) + "," + to_string(row) + "," + to_string(col) + "," + 
                to_string(toa_ns) + "," + to_string(tot_ns) + "," + to_string(cal_code) + "," + to_string(ea)
            );
        }

        if (current_event_has_ea_error) s_.events_with_ea_errors++;

        push_word_bytes(event_bytes_, pending_trailer_w_);
        bool is_crc_mismatch = (etroc::calculate_crc8(event_bytes_) != 0x00);

        s_.expected_total_hits += trl_h;
        bool has_missing_hits = (trl_h > event_hits);
        if (has_missing_hits) {
            s_.events_hit_mismatch++; s_.sw_lost_hits += (trl_h - event_hits);
        } else s_.events_perfect_hits++;

        if (!is_crc_mismatch) s_.crc_match++;
        else {
            s_.crc_mismatch++;
            bool is_silent = true;
            if (current_event_has_ea_error) { s_.crc_err_ea++; is_silent = false; }
            if (has_missing_hits) { s_.crc_err_missing++; is_silent = false; }
            if (has_missing_hits && current_event_has_ea_error) s_.crc_err_missing_with_ea++;
            if (is_silent) s_.crc_err_silent++;
        }

        if (has_missing_hits) {
            s_.events_with_missing_hits++;
            if (is_crc_mismatch) s_.events_missing_hits_and_crc_err++;
            else s_.events_missing_hits_but_crc_ok++;
        }
        if (current_event_has_ea_error) {
            s_.events_with_ea_errors++;
            if (is_crc_mismatch) s_.events_with_ea_and_crc_err++;
        }
        event_id_++;
    }

public:
    EtrocChannelAnalyzer(const string& name, uint32_t chip_id, const string& prefix, uint8_t chan_tag)
        : name_(name), expected_chip_id_(chip_id), file_prefix_(prefix), chan_tag_(chan_tag) {
        string t_name = name_ + "_" + file_prefix_;
        h_toa_ = new TH1D(Form("h_toa_%s", t_name.c_str()), Form("TOA Distribution (%s);TOA (ns);Counts", name_.c_str()), 250, 0, 25);
        h_tot_ = new TH1D(Form("h_tot_%s", t_name.c_str()), Form("TOT Distribution (%s);TOT (ns);Counts", name_.c_str()), 250, 0, 25);
        h_cal_ = new TH1D(Form("h_cal_%s", t_name.c_str()), Form("CAL Code Distribution (%s);CAL Raw Code [0-1023];Counts", name_.c_str()), 1024, 0, 1024);
    }

    ~EtrocChannelAnalyzer() { 
        delete h_toa_; delete h_tot_; delete h_cal_;
    }

    void process_word(uint64_t w, uint64_t file_total_records) {
        global_word_index_++;
        recent_words_.push_back(w);
        if (recent_words_.size() > 10) recent_words_.erase(recent_words_.begin());

        s_.total_frames++;
        FrameType payload_t = get_payload_type(w);

        if (payload_t == TRAILER && etroc::trl_chipid(w) != expected_chip_id_) {
            s_.fake_trailers_blocked++; payload_t = UNKNOWN;
        }
        if (payload_t == UNKNOWN) s_.unknown_frames++;

        if (!is_synced_) {
            if (payload_t == HEADER) is_synced_ = true;
            else { s_.words_dropped_at_start++; return; }
        }

        if (payload_t == DATA) {
            s_.hits_raw++;
            if (!in_event_) { s_.orphan_data_words++; return; }
            current_event_data_.push_back(w);
            push_word_bytes(event_bytes_, w);
            return;
        }

        if (payload_t == HEADER) {
            s_.headers++;
            uint8_t curr_l1_cnt = etroc::hdr_l1_cnt(w);
            bool l1_jump_detected = false;
            
            if (in_event_) { 
                s_.wp_boundary_cuts++; 
                expect_l1_jump_ = true; 
                
                if (global_word_index_ < (file_total_records - 100)) {
                    s_.weird_wp_boundary_cuts++;
                    if (cut_print_count_ < 5) { 
                        std::cout << "\n========================================================\n";
                        std::cout << " [\033[1;31mWeird WP Boundary Cut\033[0m] 위치: " << global_word_index_ << " / " << file_total_records << "\n";
                        std::cout << "  - 채널 : " << name_ << " | 이벤트 L1: " << (int)last_l1_cnt_ << "\n";
                        std::cout << "  - 현상 : 트레일러가 중간에 유실된 채 새 헤더가 난입했습니다.\n";
                        for (size_t i = 0; i < recent_words_.size(); i++) {
                            uint64_t rw = recent_words_[i];
                            std::string type_str = "UNKNOWN";
                            if (etroc::is_header(rw)) type_str = "HEADER ";
                            else if (etroc::is_trailer(rw)) type_str = "TRAILER";
                            else if (etroc::is_filler(rw)) type_str = "FILLER ";
                            else if (etroc::is_data(rw)) type_str = "DATA   ";
                            
                            std::cout << "    " << (i == recent_words_.size() - 1 ? ">>> " : "    ") 
                                      << "0x" << std::setfill('0') << std::setw(10) << std::hex << rw << std::dec 
                                      << " | " << type_str << "\n";
                        }
                        std::cout << "========================================================\n";
                        cut_print_count_++;
                    }
                }
            }

            if (!is_first_event_ && !expect_l1_jump_) {
                uint8_t std_cnt = (last_l1_cnt_ + 1) & 0xFF;
                if (curr_l1_cnt != std_cnt) {
                    s_.l1_counter_jumps++; s_.l1_jump_sum += ((curr_l1_cnt - last_l1_cnt_) & 0xFF);
                    l1_jump_detected = true;
                }
            }
            is_first_event_ = false; last_l1_cnt_ = curr_l1_cnt; expect_l1_jump_ = false;

            if (has_pending_event_) {
                if (l1_jump_detected) s_.events_discarded_by_jump++;
                else commit_pending_event();
                has_pending_event_ = false;
            }
            in_event_ = true; current_event_data_.clear(); event_bytes_.clear();
            push_word_bytes(event_bytes_, w); hdr_bcid_ = etroc::hdr_bcid(w); hdr_type_ = etroc::hdr_type(w);
            return;
        }

        if (payload_t == TRAILER) {
            s_.trailers++;
            if (!in_event_) { s_.events_missing_header++; return; }
            pending_trailer_w_ = w; has_pending_event_ = true; in_event_ = false; return;
        }
        if (payload_t == FILLER) s_.fillers++;
    }

    void finalize_stream(SystemStats& sys) {
        if (has_pending_event_) { commit_pending_event(); has_pending_event_ = false; }
        if (in_event_) { s_.wp_boundary_cuts++; in_event_ = false; } 
        g_stats[chan_tag_].add(s_);
        g_sys[chan_tag_].add(sys);
    }

    void export_python_data(const string& base_filename, const SystemStats& sys) const {
        string prefix = base_filename.substr(0, base_filename.find(".dat"));
        if (prefix.find(".bin") != string::npos) prefix = prefix.substr(0, prefix.find(".bin"));
        size_t last_slash = prefix.find_last_of('/');
        if (last_slash != string::npos) prefix = prefix.substr(last_slash + 1);
        string out_path_prefix = "results/" + prefix;

        // [수정됨] 개별 파일 분석 시 TCanvas를 이용한 .png 그림 생성 코드를 완전히 제거했습니다.
        // 디렉토리 오염을 방지하고 순수 대용량 고속 처리가 가능해집니다.

        ofstream hsf(out_path_prefix + "_" + name_ + "_hit_stream.csv");
        hsf << "EventID,Row,Col,TOA_ns,TOT_ns,CAL_Code,EA_Code\n";
        for(auto& str : hit_stream_buffer_) hsf << str << "\n";
        hsf.close();

        ofstream hmf(out_path_prefix + "_" + name_ + "_heatmap.csv");
        for (int r = 0; r < 16; r++) {
            for (int c = 0; c < 16; c++) hmf << s_.hit_map[r][c] << (c == 15 ? "" : ",");
            hmf << "\n";
        }

        double bufempty_ratio = sys.total_32b_words > 0 ? (double)sys.total_ram_buffer_empty / sys.total_32b_words * 100.0 : 0.0;
        double l1_jump_avg = s_.l1_counter_jumps > 0 ? (double)s_.l1_jump_sum / s_.l1_counter_jumps : 0.0;
        uint64_t total_complete = s_.events_perfect_hits + s_.events_hit_mismatch + s_.events_corrupted_bcid + s_.events_discarded_by_jump;

        double avg_toa = s_.valid_hits > 0 ? h_toa_->GetMean() : 0.0;
        double rms_toa = s_.valid_hits > 0 ? h_toa_->GetRMS() : 0.0;
        double avg_tot = s_.valid_hits > 0 ? h_tot_->GetMean() : 0.0;
        double rms_tot = s_.valid_hits > 0 ? h_tot_->GetRMS() : 0.0;
        double avg_cal = s_.valid_hits > 0 ? h_cal_->GetMean() : 0.0; 
        double rms_cal = s_.valid_hits > 0 ? h_cal_->GetRMS() : 0.0;  

        ofstream stf(out_path_prefix + "_" + name_ + "_stats.csv");
        stf << "Category,Metric,Value\n"
            << "FRAME,Total_Frames," << s_.total_frames << "\n"
            << "FRAME,Raw_Hit_Frames," << s_.hits_raw << "\n"
            << "FRAME,Headers," << s_.headers << "\n"
            << "FRAME,Trailers," << s_.trailers << "\n"
            << "FRAME,Fillers," << s_.fillers << "\n"
            << "FRAME,Fake_Trailers_Blocked," << s_.fake_trailers_blocked << "\n"
            << "FRAME,Unknown_Frames," << s_.unknown_frames << "\n"
            << "EXCLUDE,WP_Boundary_Cuts," << s_.wp_boundary_cuts << "\n"
            << "EXCLUDE,Weird_WP_Boundary_Cuts," << s_.weird_wp_boundary_cuts << "\n"
            << "EXCLUDE,Words_Dropped_At_Start," << s_.words_dropped_at_start << "\n"
            << "EXCLUDE,Events_Missing_Header," << s_.events_missing_header << "\n"
            << "EXCLUDE,Orphan_Data_Words," << s_.orphan_data_words << "\n"
            << "EVENT,Total_Complete_Events," << total_complete << "\n"
            << "EVENT,Events_Perfect_Hits," << s_.events_perfect_hits << "\n"
            << "EVENT,Events_Corrupted_BCID," << s_.events_corrupted_bcid << "\n"
            << "EVENT,Events_Hit_Mismatch," << s_.events_hit_mismatch << "\n"
            << "EVENT,Events_Discarded_By_Jump," << s_.events_discarded_by_jump << "\n"
            << "EVENT,L1_Counter_Jumps," << s_.l1_counter_jumps << "\n"
            << "EVENT,L1_Jump_Avg," << l1_jump_avg << "\n"
            << "EVENT,CRC_Match," << s_.crc_match << "\n"
            << "EVENT,CRC_Mismatch," << s_.crc_mismatch << "\n"
            << "HIT,Expected_Total_Hits," << s_.expected_total_hits << "\n"
            << "HIT,Valid_Hits," << s_.valid_hits << "\n"
            << "HIT,Lost_Hits," << s_.sw_lost_hits << "\n"
            << "HIT,Avg_TOA_ns," << avg_toa << "\n"
            << "HIT,RMS_TOA_ns," << rms_toa << "\n"
            << "HIT,Avg_TOT_ns," << avg_tot << "\n"
            << "HIT,RMS_TOT_ns," << rms_tot << "\n"
            << "HIT,Avg_CAL_Code," << avg_cal << "\n"
            << "HIT,RMS_CAL_Code," << rms_cal << "\n"
            << "HIT,EA_00_Clean," << s_.ea_counts[0] << "\n"
            << "HIT,EA_01_Corrected," << s_.ea_counts[1] << "\n"
            << "HIT,EA_10_Uncorrectable," << s_.ea_counts[2] << "\n"
            << "HIT,EA_11_Undefined," << s_.ea_counts[3] << "\n"
            << "SYSTEM,RAM_Buffer_Empty_Ratio," << bufempty_ratio << "\n";
    }
};

// =====================================================================
// 3. Raw Data Record Builder
// =====================================================================
struct RawRecord { uint64_t d; uint8_t chan_tag; uint8_t hw_type; };

vector<RawRecord> read_raw_dat(const string& path) {
    vector<RawRecord> out;
    ifstream f(path, ios::binary);
    if (!f.is_open()) return out;
    f.seekg(0, ios::end); size_t n = (size_t)f.tellg() / 8; f.seekg(0, ios::beg);
    out.reserve(n); uint8_t b[8];
    for (size_t i = 0; i < n; i++) {
        f.read((char*)b, 8); if (!f) break;
        uint64_t d = ((uint64_t)b[4]<<32)|((uint64_t)b[3]<<24)|((uint64_t)b[2]<<16)|((uint64_t)b[1]<<8)|(uint64_t)b[0];
        out.push_back({d, b[5], b[6]});
    }
    return out;
}

void run_raw_dat(const string& path) {
    vector<RawRecord> recs = read_raw_dat(path);
    cout << "Processing: " << path << "\n";

    map<int,uint64_t> tag_counts;
    for (auto& r : recs) tag_counts[r.chan_tag]++;
    vector<int> chan_tags;
    for (auto& kv : tag_counts) chan_tags.push_back(kv.first);

    size_t last_slash = path.find_last_of('/');
    string prefix = (last_slash == string::npos) ? path : path.substr(last_slash + 1);

    for (int tag : chan_tags) {
        if (tag != 0x0c && tag != 0x2a) continue;
        string label = (tag == 0x0c) ? "CH_A" : "CH_B";
        
        if (g_toa.find(tag) == g_toa.end()) {
            g_toa[tag] = new TH1D(Form("g_toa_%d", tag), Form("Global TOA (%s);TOA (ns);Counts", label.c_str()), 250, 0, 25);
            g_tot[tag] = new TH1D(Form("g_tot_%d", tag), Form("Global TOT (%s);TOT (ns);Counts", label.c_str()), 250, 0, 25);
            g_cal[tag] = new TH1D(Form("g_cal_%d", tag), Form("Global CAL Code (%s);CAL Raw Code [0-1023];Counts", label.c_str()), 1024, 0, 1024);
        }

        map<uint32_t, uint64_t> chipid_votes;
        for (auto& r : recs) {
            if (r.chan_tag != tag) continue;
            if (r.hw_type == 0x02 && (((r.d >> 9) & 0xFFFF) != 0)) chipid_votes[(r.d >> 22) & 0x1FFFF]++;
        }
        uint32_t detected_chip_id = 0x1ABCD; uint64_t best_votes = 0;
        for (auto& kv : chipid_votes) if (kv.second > best_votes) { best_votes = kv.second; detected_chip_id = kv.first; }

        EtrocChannelAnalyzer ana(label, detected_chip_id, prefix, tag);
        
        uint64_t total_recs = recs.size();
        for (auto& r : recs) {
            if (r.chan_tag != tag) continue;
            if (r.hw_type == 0x02 && ((r.d >> 9) & 0xFFFF) == 0) continue;
            ana.process_word(r.d, total_recs);
        }
        
        SystemStats sys{};
        ana.finalize_stream(sys);
        ana.export_python_data(path, sys);
    }
}

// 글로벌 통계 및 누적 플롯 생성 (results/ 폴더 고정)
void export_global_stats() {
    for (auto& kv : g_stats) {
        uint8_t tag = kv.first;
        if (tag != 0x0c && tag != 0x2a) continue;
        string label = (tag == 0x0c) ? "CH_A" : "CH_B";
        string out_path = "results/global_" + label;

        RecoStats& s = kv.second;
        SystemStats& sys = g_sys[tag];

        ofstream hmf(out_path + "_heatmap.csv");
        for (int r = 0; r < 16; r++) {
            for (int c = 0; c < 16; c++) hmf << s.hit_map[r][c] << (c == 15 ? "" : ",");
            hmf << "\n";
        }

        TCanvas* c_glob = new TCanvas("c_glob", "Global", 800, 600);
        gStyle->SetOptStat(111111);
        if (s.valid_hits > 0) {
            // 1. TOA 글로벌 플롯 & 정밀 피크 줌인 플롯 (±1.5 ns)
            g_toa[tag]->GetXaxis()->SetRangeUser(0, 25);
            g_toa[tag]->SetFillColor(38); g_toa[tag]->Draw(); c_glob->SaveAs((out_path + "_TOA_ns.png").c_str());
            
            double g_m_toa = g_toa[tag]->GetMean();
            g_toa[tag]->GetXaxis()->SetRangeUser(g_m_toa - 1.5, g_m_toa + 1.5); 
            g_toa[tag]->Draw(); c_glob->SaveAs((out_path + "_TOA_ns_zoom.png").c_str());

            // 2. TOT 글로벌 플롯 & 정밀 피크 줌인 플롯 (±1.5 ns)
            g_tot[tag]->GetXaxis()->SetRangeUser(0, 25);
            g_tot[tag]->SetFillColor(46); g_tot[tag]->Draw(); c_glob->SaveAs((out_path + "_TOT_ns.png").c_str());
            
            double g_m_tot = g_tot[tag]->GetMean();
            g_tot[tag]->SetFillColor(46);
            // [수정된 부분] g_tot_mean 대신 g_m_tot 사용!
            g_tot[tag]->GetXaxis()->SetRangeUser(g_m_tot - 1.5, g_m_tot + 1.5);
            g_tot[tag]->Draw(); c_glob->SaveAs((out_path + "_TOT_ns_zoom.png").c_str());

            // 3. CAL 글로벌 플롯 & 정밀 줌인 플롯 (±30 Bins)
            g_cal[tag]->GetXaxis()->SetRangeUser(0, 1024);
            g_cal[tag]->SetFillColor(30); g_cal[tag]->Draw(); c_glob->SaveAs((out_path + "_CAL_Code.png").c_str());
            
            double mean_cal = g_cal[tag]->GetMean();
            g_cal[tag]->GetXaxis()->SetRangeUser(mean_cal - 30.0, mean_cal + 30.0);
            g_cal[tag]->Draw(); c_glob->SaveAs((out_path + "_CAL_Code_zoom.png").c_str());
        }
        delete c_glob;

        double bufempty_ratio = sys.total_32b_words > 0 ? (double)sys.total_ram_buffer_empty / sys.total_32b_words * 100.0 : 0.0;
        double l1_jump_avg = s.l1_counter_jumps > 0 ? (double)s.l1_jump_sum / s.l1_counter_jumps : 0.0;
        uint64_t total_complete = s.events_perfect_hits + s.events_hit_mismatch + s.events_corrupted_bcid + s.events_discarded_by_jump;

        ofstream stf(out_path + "_stats.csv");
        stf << "Category,Metric,Value\n"
            << "FRAME,Total_Frames," << s.total_frames << "\n"
            << "FRAME,Raw_Hit_Frames," << s.hits_raw << "\n"
            << "FRAME,Headers," << s.headers << "\n"
            << "FRAME,Trailers," << s.trailers << "\n"
            << "FRAME,Fillers," << s.fillers << "\n"
            << "FRAME,Fake_Trailers_Blocked," << s.fake_trailers_blocked << "\n"
            << "FRAME,Unknown_Frames," << s.unknown_frames << "\n"
            << "EXCLUDE,WP_Boundary_Cuts," << s.wp_boundary_cuts << "\n"
            << "EXCLUDE,Weird_WP_Boundary_Cuts," << s.weird_wp_boundary_cuts << "\n"
            << "EXCLUDE,Words_Dropped_At_Start," << s.words_dropped_at_start << "\n"
            << "EXCLUDE,Events_Missing_Header," << s.events_missing_header << "\n"
            << "EXCLUDE,Orphan_Data_Words," << s.orphan_data_words << "\n"
            << "EVENT,Total_Complete_Events," << total_complete << "\n"
            << "EVENT,Events_Perfect_Hits," << s.events_perfect_hits << "\n"
            << "EVENT,Events_Corrupted_BCID," << s.events_corrupted_bcid << "\n"
            << "EVENT,Events_Hit_Mismatch," << s.events_hit_mismatch << "\n"
            << "EVENT,Events_Discarded_By_Jump," << s.events_discarded_by_jump << "\n"
            << "EVENT,L1_Counter_Jumps," << s.l1_counter_jumps << "\n"
            << "EVENT,L1_Jump_Avg," << l1_jump_avg << "\n"
            << "EVENT,CRC_Match," << s.crc_match << "\n"
            << "EVENT,CRC_Mismatch," << s.crc_mismatch << "\n"
            << "HIT,Expected_Total_Hits," << s.expected_total_hits << "\n"
            << "HIT,Valid_Hits," << s.valid_hits << "\n"
            << "HIT,Lost_Hits," << s.sw_lost_hits << "\n"
            << "HIT,Avg_TOA_ns," << (s.valid_hits > 0 ? g_toa[tag]->GetMean() : 0.0) << "\n"
            << "HIT,RMS_TOA_ns," << (s.valid_hits > 0 ? g_toa[tag]->GetRMS() : 0.0) << "\n"
            << "HIT,Avg_TOT_ns," << (s.valid_hits > 0 ? g_tot[tag]->GetMean() : 0.0) << "\n"
            << "HIT,RMS_TOT_ns," << (s.valid_hits > 0 ? g_tot[tag]->GetRMS() : 0.0) << "\n"
            << "HIT,Avg_CAL_Code," << (s.valid_hits > 0 ? g_cal[tag]->GetMean() : 0.0) << "\n"
            << "HIT,RMS_CAL_Code," << (s.valid_hits > 0 ? g_cal[tag]->GetRMS() : 0.0) << "\n"
            << "HIT,EA_00_Clean," << s.ea_counts[0] << "\n"
            << "HIT,EA_01_Corrected," << s.ea_counts[1] << "\n"
            << "HIT,EA_10_Uncorrectable," << s.ea_counts[2] << "\n"
            << "HIT,EA_11_Undefined," << s.ea_counts[3] << "\n"
            << "SYSTEM,RAM_Buffer_Empty_Ratio," << bufempty_ratio << "\n";
    }
}

int main(int argc, char** argv) {
    gROOT->SetBatch(kTRUE);
    if (argc < 2) { cout << "Usage: " << argv[0] << " raw_data/output_run_*_rb0.dat\n"; return 1; }

    cout << "\n========================================================\n";
    cout << "  Batch parsing mode active. Outputs routing to results/\n";
    cout << "========================================================\n";

    for (int i = 1; i < argc; i++) run_raw_dat(argv[i]);
    
    cout << "\nExporting Global Summaries...\n";
    
    ofstream summary_csv("results/global_timing_summary.csv");
    summary_csv << "Channel,Total_Valid_Hits,Avg_TOA_ns,RMS_TOA_ns,Avg_TOT_ns,RMS_TOT_ns,Avg_CAL_Code,RMS_CAL_Code\n";
    
    for (auto& kv : g_toa) {
        uint8_t tag = kv.first;
        if (tag != 0x0c && tag != 0x2a) continue;
        string label = (tag == 0x0c) ? "CH_A" : "CH_B";
        
        double mean_toa = g_valid_hits[tag] > 0 ? g_toa[tag]->GetMean() : 0;
        double rms_toa  = g_valid_hits[tag] > 0 ? g_toa[tag]->GetRMS() : 0;
        double mean_tot = g_valid_hits[tag] > 0 ? g_tot[tag]->GetMean() : 0;
        double rms_tot  = g_valid_hits[tag] > 0 ? g_tot[tag]->GetRMS() : 0;
        double mean_cal = g_valid_hits[tag] > 0 ? g_cal[tag]->GetMean() : 0;
        double rms_cal  = g_valid_hits[tag] > 0 ? g_cal[tag]->GetRMS() : 0;

        summary_csv << label << "," << g_valid_hits[tag] << ","
                    << fixed << setprecision(3) << mean_toa << "," << rms_toa << ","
                    << mean_tot << "," << rms_tot << ","
                    << mean_cal << "," << rms_cal << "\n";
    }
    summary_csv.close();

    export_global_stats();
    cout << "Done! All results safely landed in the 'results/' folder.\n\n";
    
    return 0;
}
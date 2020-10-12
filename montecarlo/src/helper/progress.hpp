#ifndef _progress_h_
#define _progress_h_

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

// helper class to monitor progress of the loops
class progress_bar {
 private:
  const int   _barWidth = 100;
  int         _pos = -1;
  float       _progress = 0;
  std::string _title = "";
  int         _i_max = 0;
  int         _i = 0;

  bool _is_silent = false;

  bool        _start_chrono = true;
  std::time_t _start_time;
  std::time_t _current_time;

  const char full = '=';
  const char empty = '.';

 public:
  // constructor to initialize internal state of the progress bar
  progress_bar(const int& i_max, const std::string& title = "", const bool& is_silent = false) {
    _title = title;
    _i_max = i_max - 1;
    _i = 0;

    _is_silent = is_silent;

    if (not _is_silent) {
      std::cout << _title << ":" << std::endl;
    }
  };

  std::string estimate_remaining_time() {
    std::string remaining_time;
    if (_start_chrono) {
      _start_time = std::time(nullptr);
      _current_time = std::time(nullptr);
      _start_chrono = false;
    } else {
      _current_time = std::time(nullptr);
    }

    int sec = int(std::difftime(_current_time, _start_time) * (1 - _progress) / (_progress + 0.001));

    int hour = sec / 3600;
    sec = sec % 3600;
    int min = sec / 60;
    sec = sec % 60;
    std::stringstream ss;
    ss << "remaining time " << std::setw(2) << std::setfill('0') << hour << ":" << std::setw(2) << std::setfill('0')
       << min << ":" << std::setw(2) << std::setfill('0') << sec;
    return ss.str();
  };

  // stepping function to keep track of progress internally
  void step() {
    if (_is_silent) return;

    _progress = float(_i) / float(_i_max);
    _i++;
    

    int new_pos = _progress * _barWidth;

    if (new_pos == _pos) return;

    _pos = new_pos;
    std::cout << "[";
    for (int j = 0; j < _barWidth; ++j) {
      if (j < _pos)
        std::cout << full;
      else
        std::cout << empty;
    }
    std::cout << "] " << int((_progress)*100.0) << "% " << estimate_remaining_time() << "\r" << std::flush;
    if (_i >= _i_max) {
      std::cout << std::endl;
    }
  };

  // stepping function to pass progress level explicitly and do the rest internally
  void step(const int& i) {
    if (_is_silent) return;

    _progress = float(i) / float(_i_max);
    _pos = _progress * _barWidth;
    std::cout << "[";
    for (int j = 0; j < _barWidth; ++j) {
      if (j <= _pos)
        std::cout << full;
      else
        std::cout << empty;
    }
    std::cout << "] " << int((_progress)*100.0) << "% " << estimate_remaining_time() << "\r" << std::flush;
    if (i == _i_max) {
      std::cout << std::endl;
    }
  };

  // constructor without any internal state initialization
  progress_bar(const bool& is_silent = false) { _is_silent = is_silent; };

  // stepping function that does not need initialization
  void step(const int& i, const int& i_max, const std::string& title) {
    if (_is_silent) return;

    _progress = float(i) / float(i_max - 1);
    _pos = _progress * _barWidth;
    std::cout << title << ": [";
    for (int j = 0; j < _barWidth; ++j) {
      if (j <= _pos)
        std::cout << full;
      else
        std::cout << empty;
    }
    std::cout << "] " << int((_progress)*100.0) << "% " << estimate_remaining_time() << "\r" << std::flush;
    if (i == i_max - 1) {
      std::cout << std::endl;
    }
  };
};

#endif  //_progress_h_
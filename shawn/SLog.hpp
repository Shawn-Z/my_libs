#ifndef SHAWN_LOG_HPP
#define SHAWN_LOG_HPP

#include <vector>
#include <ctime>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <glog/logging.h>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h>

namespace shawn {
    class SLog {
    public:
        void init(const char* p_name, std::string p_folder = "log_dir", google::LogSeverity p_stderr_severity = google::ERROR) {
            google::InitGoogleLogging(p_name);
            std::string folder = getenv("HOME");
            folder += "/";
            folder += p_folder;
            if (access(folder.c_str(), 06) < 0) {
                while (mkdir(folder.c_str(), 0777) < 0) {
                    std::cerr<<"glog -- create folder failed"<<std::endl;
                }
            }
            FLAGS_log_dir = folder;
            google::SetStderrLogging(p_stderr_severity);
        }
    };
}

#endif //SHAWN_LOG_HPP
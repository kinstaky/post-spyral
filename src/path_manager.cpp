#include "include/path_manager.h"

#include <toml.hpp>

namespace gong {

PathManager& PathManager::Instance() {
	static PathManager instance;
	return instance;
}

void PathManager::Initialize(const char *path) {
	std::lock_guard<std::mutex> lock(mutex_);
	exec_path_ = (std::filesystem::current_path() / path).parent_path();
	// read config.toml
	std::filesystem::path config_file_path = exec_path_ / "../../config.toml";
	// read config.toml
	auto config = toml::parse(config_file_path);
	// get data path
	data_path_ = config["path"]["data"].as_string();
	// make sure data path is absolute
	data_path_ = std::filesystem::absolute(data_path_);

	init_ = true;
}

std::filesystem::path PathManager::ExecPath() const {
	std::lock_guard<std::mutex> lock(mutex_);
	if (!init_) throw std::runtime_error("PathManager::ExecPath() called before Initialize()");
	return exec_path_.lexically_normal();
}


std::filesystem::path PathManager::DataPath() const {
	std::lock_guard<std::mutex> lock(mutex_);
	if (!init_) throw std::runtime_error("PathManager::DataPath() called before Initialize()");
	return data_path_.lexically_normal();
}


}	// namespace gong
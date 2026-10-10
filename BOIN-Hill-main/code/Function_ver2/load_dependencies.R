load_all_r_files <- function(directory) {
  # 디렉토리 내의 모든 .R 파일 검색
  r_files <- list.files(
    path = directory,
    pattern = "\\.R$",
    full.names = TRUE,
    recursive = TRUE  # 하위 디렉토리도 포함
  )
  
  # 로드된 파일 추적
  loaded_files <- character(0)
  
  # 소스 파일 내용에서 source() 호출 찾기
  find_source_calls <- function(file_path) {
    # 파일 내용 읽기
    file_content <- readLines(file_path, warn = FALSE)
    
    # source() 호출 패턴 찾기
    source_pattern <- "source\\([\"'](.+?)[\"'].*\\)"
    matches <- gregexpr(source_pattern, file_content)
    
    # 찾은 source 호출에서 파일 경로 추출
    source_files <- character(0)
    for (i in seq_along(file_content)) {
      match <- regmatches(file_content[i], gregexpr(source_pattern, file_content[i]))
      if (length(match[[1]]) > 0) {
        # 파일 경로 부분만 추출
        for (m in match[[1]]) {
          file_ref <- sub(source_pattern, "\\1", m)
          source_files <- c(source_files, file_ref)
        }
      }
    }
    
    return(source_files)
  }
  
  # 재귀적으로 파일 로드하기
  load_file <- function(file_path, parent_env = parent.frame()) {
    # 이미 로드된 파일은 건너뛰기
    if (file_path %in% loaded_files) {
      message(paste("Skip already loaded:", file_path))
      return()
    }
    
    # 파일이 존재하는지 확인
    if (!file.exists(file_path)) {
      # 상대 경로 처리
      possible_paths <- c(
        file.path(directory, file_path),
        file.path(getwd(), file_path)
      )
      
      found <- FALSE
      for (path in possible_paths) {
        if (file.exists(path)) {
          file_path <- path
          found <- TRUE
          break
        }
      }
      
      if (!found) {
        warning(paste("File not found:", file_path))
        return()
      }
    }
    
    # 로드된 파일 목록에 추가
    loaded_files <<- c(loaded_files, file_path)
    
    # 먼저 이 파일에서 source()로 참조하는 다른 파일들 찾기
    referenced_files <- find_source_calls(file_path)
    
    # 참조된 파일들 먼저 로드
    for (ref_file in referenced_files) {
      load_file(ref_file, parent_env)
    }
    
    # 이제 현재 파일 로드
    message(paste("Loading:", file_path))
    source(file_path, local = parent_env)
  }
  
  # 모든 파일 로드
  for (file in r_files) {
    load_file(file)
  }
  
  # 로드된 파일 수 반환
  return(length(unique(loaded_files)))
}

# 병렬 처리를 위한 필수 함수들을 로드하는 함수
load_essential_functions <- function() {
  # Function_ver2 디렉토리의 모든 R 파일 로드
  num_files <- load_all_r_files("./Function_ver2")
  message(paste("Total", num_files, "R files loaded."))
}
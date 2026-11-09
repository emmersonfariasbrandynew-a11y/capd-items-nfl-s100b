# EXECUTAR_ANALISE.R
# Runner: verifies the SHA-256 checksums in MANIFESTO_EXECUCAO.csv, then runs scripts 00-09 in order
# with per-step checkpoints and resume.
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

Sys.setlocale('LC_CTYPE','English_United States.utf8')
lib<-Sys.getenv('CAPD_R_LIB','')
if(nzchar(lib)) .libPaths(c(lib,.libPaths()))
stopifnot(requireNamespace('digest',quietly=TRUE))
m<-read.csv('MANIFESTO_EXECUCAO.csv',stringsAsFactors=FALSE)
for(i in seq_len(nrow(m)))stopifnot(file.exists(m$arquivo[i]),identical(digest::digest(file=m$arquivo[i],algo='sha256'),m$sha256[i]))
manifest_hash<-digest::digest(file='MANIFESTO_EXECUCAO.csv',algo='sha256')
run<-new.env(parent=globalenv())
fs<-c('00_setup.R','00b_motor_deltarho.R','01_data_prep.R','02_descritivas.R','03_precisao_MDE.R','04_S1_primaria.R','05_S2_dominios.R','06_H2H3_diferenca_rho.R','07_H5_comparadores.R','08_transicoes.R','09_H4_sensibilidades.R')
dir.create('output',showWarnings=FALSE)
status_path<-'output/STATUS_EXECUCAO.csv'
st<-data.frame(script=character(),estado=character(),inicio=character(),fim=character(),erro=character())
if(file.exists(status_path))st<-read.csv(status_path,stringsAsFactors=FALSE,na.strings='')
completed<-character()
if(file.exists('output/CHECKPOINT_ETAPA.rds')){
 cp<-readRDS('output/CHECKPOINT_ETAPA.rds');stopifnot(identical(cp$manifest_hash,manifest_hash))
 list2env(cp$objects,envir=run);completed<-cp$completed
 if(!is.null(cp$rng))assign('.Random.seed',cp$rng,envir=.GlobalEnv)
}
capture.output(sessionInfo(),file='output/SESSION_INFO_INICIO.txt')
for(s in fs){
 if(s %in% completed){cat('Etapa concluída recuperada:',s,'\n');next}
 st<-rbind(st,data.frame(script=s,estado='EM_EXECUCAO',inicio=as.character(Sys.time()),fim='',erro=''))
 write.csv(st,status_path,row.names=FALSE)
 cat('\nINICIANDO',s,as.character(Sys.time()),'\n');flush.console()
 err<-tryCatch({source(file.path('scripts',s),local=run,encoding='UTF-8');NULL},error=function(e)conditionMessage(e))
 if(!is.null(err)){
  st$estado[nrow(st)]<-'FALHOU';st$erro[nrow(st)]<-err;st$fim[nrow(st)]<-as.character(Sys.time());write.csv(st,status_path,row.names=FALSE)
  stop(err)
 }
 if(exists('CFG',envir=run,inherits=FALSE))stopifnot(run$CFG$n_boot==5000,run$CFG$seed==20260818)
 completed<-c(completed,s)
 objs<-as.list(run,all.names=TRUE);objs<-objs[!vapply(objs,inherits,logical(1),'connection')]
 saveRDS(list(manifest_hash=manifest_hash,completed=completed,objects=objs,rng=if(exists('.Random.seed',envir=.GlobalEnv))get('.Random.seed',envir=.GlobalEnv) else NULL),'output/CHECKPOINT_ETAPA.tmp.rds')
 if(!file.copy('output/CHECKPOINT_ETAPA.tmp.rds','output/CHECKPOINT_ETAPA.rds',overwrite=TRUE))stop('Falha ao preservar checkpoint')
 st$estado[nrow(st)]<-'CONCLUIDO';st$fim[nrow(st)]<-as.character(Sys.time());write.csv(st,status_path,row.names=FALSE)
 cat('CONCLUIDO',s,as.character(Sys.time()),'\n');flush.console()
}
capture.output(sessionInfo(),file='output/SESSION_INFO_FINAL.txt')
writeLines(c('EXECUCAO FINAL CONCLUIDA',as.character(Sys.time()),'Revisar avisos e contratos de saída antes de interpretar.'),'output/CONCLUSAO_EXECUCAO.txt')

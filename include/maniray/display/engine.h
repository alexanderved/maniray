#ifndef _MR_ENGINE_H
#define _MR_ENGINE_H

#include "window.h"
#include "program.h"

typedef struct mr_engine mr_engine;

typedef void (*mr_engine_update)(mr_engine *, void *);

int mr_engine_init(mr_window *window);

mr_engine *mr_engine_create(mr_window *window);
void mr_engine_destroy(mr_engine *engine);

void mr_engine_run(mr_engine *engine, mr_program *compute_program, mr_engine_update update, void *userdata);

#endif // _MR_ENGINE_H
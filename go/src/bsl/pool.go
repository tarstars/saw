package bsl

import "sync"

//Struct for pool of processes
type Pool struct {
	taskChannel chan Runnable
	waitGroup sync.WaitGroup
}

//Thread of executing tasks
func (pool *Pool) executor() {
	for {
		if task, ok := <- pool.taskChannel; ok {
			task.Run()
			pool.waitGroup.Done()
		} else {
			return
		}
	}
}

//Interface for tasks
type Runnable interface {
	Run()
}

//Create new pool for n async processes
func NewPool(n int) (ret *Pool) {
	ret = &Pool{}
	ret.taskChannel = make(chan Runnable)
	for t:=0; t < n; t++ {
		go ret.executor()
	}
	return
}

//Add task in pool
func (pool *Pool) AddTask(task Runnable) {
	pool.waitGroup.Add(1)
	pool.taskChannel <- task
}

//Wait for all task finish
func (pool *Pool) WaitAll() {
	pool.waitGroup.Wait()
}

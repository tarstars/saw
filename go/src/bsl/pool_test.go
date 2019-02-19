package bsl

import (
	"testing"
)

//task for set zero in given index
type SetOne struct {
	a []int
	ind int
}

//set zero in array
func (f *SetOne) Run() {
	f.a[f.ind] = 0
}

//fill array with 1000 "1" and then change them to "0"
//in threads governed by the pool
func TestPoolTask(t *testing.T) {
	n := 1000
	a := make([]int, n)
	for t := range a {
		a[t] = 1
	}

	pool := NewPool(10)
	for t := range a {
		pool.AddTask(&SetOne{a, t})
	}
	pool.WaitAll()

	for ind, v := range a {
		if v != 0 {
			t.Error("not 0 at index", ind)
			return
		}
	}
}

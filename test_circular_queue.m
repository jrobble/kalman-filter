queue = circular_queue(3);

[queue,val] = queue.mean();
val

[queue] = queue.push_back([1,1]);
[queue,val] = queue.mean();
val

[queue] = queue.push_back([3,3]);
[queue,val] = queue.mean();
val

[queue] = queue.push_back([2,2]);
[queue,val] = queue.mean();
val

[queue] = queue.push_back([6,6]);
[queue,val] = queue.mean();
val

[queue] = queue.push_back([100,0]);
[queue,val] = queue.mean();
val
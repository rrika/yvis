# yvis: a vvis-compatible visibility pre-computer

The launcher is written in python and is responsible for putting the resulting data back into the .bsp. The main algorithm is written in Rust.

```
git clone --recursive https://github.com/rrika/yvis.git # download yvis and srctools
python yvis.py input.bsp                                # assumes existence of input.prt next to input.bsp
cargo run --release input.prt output.pvs output.ppm     # output visibility blob and ppm image for debugging
```

I hope to improve upon existing visibility tools by:
- using more accurate linear-programming based visibility checks
- using sparse data structures for performance
- using coarse-to-fine approximations to reduce work
- using performance statistics to auto tune the process
- giving hints for visleaf merging/splitting based on estimated triangles seen

For more information see my [notes of visibility](https://github.com/rrika/notes/wiki/Visibility).

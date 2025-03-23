# Compact-PIM-based-Design
This repository contains the code associated with "Optimizing and Exploring System Performance in Compact Processing-in-Memory-based Chips", accepted to AICAS2025. (https://arxiv.org/abs/2502.21259)

## Introduction
Processing-in-memory (PIM) is a promising computing paradigm to tackle the “memory wall” challenge. However, PIM system-level benefits over traditional von Neumann architecture can be reduced when the memory array cannot fully store all the neural network (NN) weights. The NN size is increasing while the PIM design size cannot scale up accordingly due to area constraints. Therefore, this work targets the system performance optimization and exploration for compact PIM designs. We first analyze the impact of data movement on compact designs. Then, we propose a novel pipeline method that maximizes the reuse of
NN weights to improve the throughput and energy efficiency of inference in compact chips. To further boost throughput, we introduce a scheduling algorithm to mitigate the pipeline bubble problem. Moreover, we investigate the trade-off between the network size and system performance for a compact PIM chip. Experimental results show that the proposed algorithm achieves 2.35× and 0.5% improvement in throughput and energy efficiency, respectively. Compared to the area-unlimited design, our compact chip achieves approximately 56.5% of the throughput and 58.6% of the energy efficiency while using only onethird of the chip area, along with 1.3× improvement in area efficiency. Our compact design also outperforms the modern GPU with 4.56× higher throughput and 157× better energy efficiency. Besides, our compact design uses less than 20% of the system energy for data movement as batch size scales up.

## Citation
```
@article{chen2025optimizing,
  title={Optimizing and Exploring System Performance in Compact Processing-in-Memory-based Chips},
  author={Chen, Peilin and Yang, Xiaoxuan},
  journal={arXiv preprint arXiv:2502.21259},
  year={2025}
}
```

## Reference Repositories
NeuroSim: https://github.com/neurosim/DNN_NeuroSim_V1.4

DRAMPower: https://github.com/tukl-msd/DRAMPower 
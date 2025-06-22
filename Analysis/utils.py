import threading
import time

def timed_input(prompt, timeout=5):
    user_input = [None]

    def get_input():
        user_input[0] = input(prompt)

    thread = threading.Thread(target=get_input)
    thread.daemon = True
    thread.start()
    thread.join(timeout)

    if thread.is_alive():
        print("\nTime expired!")
        return None
    else:
        return user_input[0]

def repeat_and_time(func, n=1, *args, **kwargs):
    durations = []

    for i in range(1, n + 1):
        print(f"\nRun {i}/{n}...")
        start_time = time.time()

        # Run the target function
        func(*args, **kwargs)

        end_time = time.time()
        elapsed = end_time - start_time
        durations.append(elapsed)

        mins, secs = divmod(elapsed, 60)
        print(f"Duration: {int(mins)} min {secs:.2f} sec")

    return durations
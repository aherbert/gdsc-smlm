/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.ij.plugins;

import java.util.ArrayList;
import java.util.Arrays;

import gdsc.smlm.utils.Pair;

/**
 * Allow processing work in stages, repeating only the stages necessary to render new results given changes to settings.
 * This class is designed to be used to allow live display of results upon settings changes by running the analysis on
 * worker threads.
 *
 * @author Alex Herbert
 */
public class Workflow<S, R>
{
	/**
	 * Default delay (in milliseconds) to use for dialog previews
	 */
	public static final long DELAY = 500;

	private class Work
	{
		public final long timeout;
		public final Pair<S, R> work;

		Work(long time, Pair<S, R> work)
		{
			if (work.a == null)
				throw new NullPointerException("Settings cannot be null");
			this.timeout = time;
			this.work = work;
		}

		Work(Pair<S, R> work)
		{
			this(0, work);
		}

		S getSettings()
		{
			return work.a;
		}

		R getResults()
		{
			return work.b;
		}
	}

	/**
	 * Allow work to be added to a FIFO stack in a synchronised manner
	 */
	private class WorkStack
	{
		// We only support a stack size of 1
		private Work work = null;

		synchronized void setWork(Work work)
		{
			this.work = work;
		}

		synchronized void addWork(Work work)
		{
			this.work = work;
			this.notify();
		}

		@SuppressWarnings("unused")
		synchronized void close()
		{
			this.work = null;
			this.notify();
		}

		synchronized Work getWork()
		{
			Work work = this.work;
			this.work = null;
			return work;
		}

		boolean isEmpty()
		{
			return work == null;
		}
	}

	public class RunnableWorker implements Runnable
	{
		private final WorkflowWorker<S, R> worker;
		private boolean running = true;
		private Work lastWork = null;
		private Work result;
		private WorkStack inbox;
		private Object[] outbox;

		RunnableWorker(WorkflowWorker<S, R> worker)
		{
			this.worker = worker;
		}

		@Override
		@SuppressWarnings("unchecked")
		public void run()
		{
			// Note: We check the condition for loop termination within the loop
			while (true)
			{
				try
				{
					Work work = null;
					synchronized (inbox)
					{
						if (inbox.isEmpty())
						{
							debug("Inbox empty, waiting ...");
							inbox.wait();
						}
						work = inbox.getWork();
						if (work != null)
							debug(" Found work");
					}
					if (work == null)
					{
						debug(" No work, stopping");
						break;
					}

					// Delay processing the work. Allows the work to be updated before we process it.
					if (work.timeout != 0)
					{
						debug(" Checking delay");
						long timeout = work.timeout;
						while (System.currentTimeMillis() < timeout)
						{
							debug(" Delaying");
							Thread.sleep(50);
							// Assume new work can be added to the inbox. Here we are peaking at the inbox
							// so we do not take ownership with synchronized
							if (inbox.work != null)
								timeout = inbox.work.timeout;
						}
						// If we intend to modify the inbox then we should take ownership
						synchronized (inbox)
						{
							if (!inbox.isEmpty())
							{
								work = inbox.getWork();
								debug(" Found updated work");
							}
						}
					}

					if (!equals(work, lastWork))
					{
						// Create a new result
						debug(" Creating new result");
						Pair<S, R> results = worker.doWork(work.work);
						result = new Work(results);
					}
					else
					{
						// Pass through the new settings with the existing results.
						// This allows a worker to ignore settings changes that do not effect
						// its results but pass the new settings to downstream workers.
						debug(" Updating existing result");
						result = new Work(new Pair<>(work.getSettings(), result.getResults()));
					}
					lastWork = work;
					// Add the result to the output
					if (outbox != null)
					{
						debug(" Posting result");
						for (int i = outbox.length; i-- > 0;)
							((WorkStack) outbox[i]).addWork(result);
					}
				}
				catch (InterruptedException e)
				{
					debug(" Interrupted, stopping");
					break;
				}

				if (!running)
				{
					debug(" Shutdown");
					break;
				}
			}
		}

		private void debug(String msg)
		{
			if (debug)
				System.out.println(worker.getClass().getSimpleName() + msg);
		}

		private boolean equals(Work work, Work lastWork)
		{
			if (lastWork == null)
				return false;

			// We must selectively compare this as not all settings changes matter.
			if (compareNulls(work.getSettings(), lastWork.getSettings()))
				return false;
			if (!worker.equalSettings(work.getSettings(), lastWork.getSettings()))
				return false;

			// We can compare these here using object references.
			// Any new results passed in will trigger equals to fail.
			boolean result = worker.equalResults(work.getResults(), lastWork.getResults());
			if (!result)
				worker.newResults();
			return result;
		}

		private boolean compareNulls(Object o1, Object o2)
		{
			if (o1 == null)
				return o2 != null;
			return o2 == null;
		}
	}

	private WorkStack inputStack = new WorkStack();
	private ArrayList<Thread> threads;
	private ArrayList<RunnableWorker> workers = new ArrayList<>();
	private long delay = 0;

	/** The debug flag. Set to true to allow print statements during operation. */
	public boolean debug = false;

	/**
	 * Adds the worker. Connect the inbox to the previous worker outbox, or the primary input if the previous is null.
	 *
	 * @param worker
	 *            the worker
	 * @return the worker id
	 */
	public int add(WorkflowWorker<S, R> worker)
	{
		return add(worker, -1);
	}

	/**
	 * Adds the worker. Connect the inbox to the previous worker outbox, or the primary input if the previous is null.
	 * <p>
	 * Use this method to add workers that can operate in parallel on the output from a previous worker.
	 *
	 * @param worker
	 *            the worker
	 * @param previous
	 *            the previous worker id
	 * @return the worker id
	 */
	public int add(WorkflowWorker<S, R> worker, int previous)
	{
		if (previous <= 0 || previous > workers.size())
			return addToChain(worker);
		else
			return addToChain(worker, workers.get(previous - 1));
	}

	/**
	 * Adds the worker. Connect the inbox to the previous worker outbox, or the primary input.
	 *
	 * @param worker
	 *            the worker
	 * @return the worker id
	 */
	private int addToChain(WorkflowWorker<S, R> worker)
	{
		if (workers.isEmpty())
		{
			return addToChain(worker, null);
		}
		else
		{
			// Chain together
			RunnableWorker previous = workers.get(workers.size() - 1);
			return addToChain(worker, previous);
		}

	}

	/**
	 * Adds the worker. Connect the inbox to the previous worker outbox, or the primary input.
	 *
	 * @param worker
	 *            the worker
	 * @param previous
	 *            the previous worker from which to take work
	 * @return the worker id
	 */
	@SuppressWarnings("unchecked")
	private int addToChain(WorkflowWorker<S, R> inputWorker, RunnableWorker previous)
	{
		RunnableWorker worker = new RunnableWorker(inputWorker);
		if (previous == null)
		{
			// Take the primary input
			worker.inbox = inputStack;
		}
		else
		{
			// Chain together
			int size = (previous.outbox == null) ? 0 : previous.outbox.length;
			if (size == 0)
				previous.outbox = new Object[1];
			else
				previous.outbox = Arrays.copyOf(previous.outbox, size + 1);
			previous.outbox[size] = new WorkStack();
			worker.inbox = (WorkStack) previous.outbox[size];
		}
		workers.add(worker);
		return workers.size();
	}

	public synchronized void start()
	{
		shutdown(true);
		threads = startWorkers(workers);
	}

	public synchronized void shutdown(boolean now)
	{
		if (threads != null)
		{
			finishWorkers(workers, threads, now);
			threads = null;
		}
	}

	public boolean isRunning()
	{
		return threads != null;
	}

	private ArrayList<Thread> startWorkers(ArrayList<RunnableWorker> workers)
	{
		ArrayList<Thread> threads = new ArrayList<>();
		for (RunnableWorker w : workers)
		{
			Thread t = new Thread(w);
			w.running = true;
			t.setDaemon(true);
			t.start();
			threads.add(t);
		}
		return threads;
	}

	private void finishWorkers(ArrayList<RunnableWorker> workers, ArrayList<Thread> threads, boolean now)
	{
		// Finish work
		for (int i = 0; i < threads.size(); i++)
		{
			Thread t = threads.get(i);
			RunnableWorker w = workers.get(i);

			if (now)
			{
				// Stop immediately any running worker
				try
				{
					t.interrupt();
				}
				catch (SecurityException e)
				{
					// We should have permission to interrupt this thread.
					e.printStackTrace();
				}
			}
			else
			{
				// Stop after the current work in the inbox
				w.running = false;

				// Notify a workers waiting on the inbox.
				// Q. How to check if the worker is sleeping?
				synchronized (w.inbox)
				{
					w.inbox.notify();
				}

				// Leave to finish their current work
				try
				{
					t.join(0);
				}
				catch (InterruptedException e)
				{
				}
			}
		}
	}

	/**
	 * Add the work settings into the workflow queue and run.
	 *
	 * @param settings
	 *            the settings
	 */
	public void run(S settings)
	{
		run(settings, null);
	}

	/**
	 * Add the work settings into the workflow queue and run.
	 *
	 * @param settings
	 *            the settings
	 * @param results
	 *            the results
	 */
	public void run(S settings, R results)
	{
		inputStack.addWork(new Work(getTimeout(), new Pair<>(settings, results)));
	}

	/**
	 * Stage the work settings into the workflow queue but do not run.
	 *
	 * @param settings
	 *            the settings
	 */
	public void stage(S settings)
	{
		stage(settings, null);
	}

	/**
	 * Stage the work settings into the workflow queue but do not run.
	 *
	 * @param settings
	 *            the settings
	 * @param results
	 *            the results
	 */
	public void stage(S settings, R results)
	{
		inputStack.setWork(new Work(getTimeout(), new Pair<>(settings, results)));
	}

	/**
	 * Checks if work is staged.
	 *
	 * @return true, if is staged
	 */
	public boolean isStaged()
	{
		return !inputStack.isEmpty();
	}

	/**
	 * Run the staged work.
	 */
	public void runStaged()
	{
		inputStack.addWork(inputStack.work);
	}

	/**
	 * Gets the delay to allow before running new work.
	 *
	 * @return the delay
	 */
	public long getDelay()
	{
		return delay;
	}

	/**
	 * Sets the delay to allow before running new work. All work is added to the queue with a timeout equals to the
	 * current system time plus the delay. If additional work is added before the delay elapses then the
	 * preceding work is ignored and the timeout reset.
	 *
	 * @param delay
	 *            the new delay
	 */
	public void setDelay(long delay)
	{
		this.delay = Math.max(0, delay);
	}

	private long getTimeout()
	{
		if (delay == 0)
			return 0;
		return System.currentTimeMillis() + delay;
	}

	/**
	 * Start preview mode. This sets the delay to the default delay time (see {@link #DELAY}). It should be called when
	 * all new work should respect the timeout delay.
	 */
	public void startPreview()
	{
		setDelay(DELAY);
	}

	/**
	 * Stop preview. This sets the delay to zero.
	 */
	public void stopPreview()
	{
		setDelay(0);
	}
}
